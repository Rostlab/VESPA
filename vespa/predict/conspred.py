#!/usr/bin/env python
"""
 Copyright (C) 2021 Rostlab
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# -*- coding: utf-8 -*-

"""
File: conspred
Conservation Prediction for embeddings
"""

# Default Imports
from collections import namedtuple
import enum
from pathlib import Path

# Lib Imports
import torch
import torch.utils.data
import torch.nn as nn
from torch.nn.utils.rnn import pad_sequence
import h5py
from tqdm import tqdm

# Module Imports
from vespa.predict.config import CACHE_DIR, DEVICE, VERBOSE


ConservationOut = namedtuple("ConservationOut", ["probability_out", "class_out"])

# Convolutional neural network (two convolutional layers)
class ConsCNN(nn.Module):
    def __init__(self):
        super(ConsCNN, self).__init__()
        n_features = 1024
        bottleneck_dim = 32
        n_classes = 9
        dropout_rate = 0.25
        self.classifier = nn.Sequential(
            nn.Conv2d(
                n_features, bottleneck_dim, kernel_size=(7, 1), padding=(3, 0)
            ),  # 7x32
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Conv2d(bottleneck_dim, n_classes, kernel_size=(7, 1), padding=(3, 0)),
        )
        self.softmax = nn.Softmax(dim=1)

    def forward(self, x):
        """
        L = protein length
        B = batch-size
        F = number of features (1024 for embeddings)
        N = number of classes (9 for conservation)
        """
        # IN: X = (B x L x F); OUT: (B x F x L, 1)
        x = x.permute(0, 2, 1).unsqueeze(dim=-1)
        Yhat_consurf = self.classifier(x)  # OUT: Yhat_consurf = (B x N x L x 1)
        # IN: (B x N x L x 1); OUT: ( B x L x N )
        Yhat_consurf = Yhat_consurf.squeeze(dim=-1)
        return Yhat_consurf

    def extract_probabilities(self, Yhat):
        return self.softmax(Yhat)

    def extract_conservation_score(self, Yhat):
        return torch.max(Yhat.data, dim=1)[1]  # get index of most reliable class prediction


class BatchCollator(object):
    def __call__(self, batch, ignore_idx=-100):
        # batch is a list of the samples returned by your __get_item__ method in your CustomDataset
        pdb_ids, X, seq_lens = zip(*batch)
        X = pad_sequence(X, batch_first=True, padding_value=ignore_idx)
        return (list(pdb_ids), X, seq_lens)


class SequenceDataset(torch.utils.data.Dataset):
    def __init__(self, h5file):
        self.interface = h5file
        self.keys = [key for key in self.interface]
        self.data_len = len(self.keys)  # number of samples in the set

    def __len__(self):
        return self.data_len

    def __getitem__(self, index):
        pdb_id = self.keys[index]
        x = torch.from_numpy(self.interface[pdb_id][:]).to(torch.float32)
        return (pdb_id, x, x.shape[0])



def get_dataloader(customdata, batch_size):
    # Create dataloaders with collate function
    collator = BatchCollator()
    dataset = SequenceDataset(customdata)
    return torch.utils.data.DataLoader(
        dataset=dataset,
        batch_size=batch_size,
        shuffle=False,
        drop_last=False,
        num_workers=3,
        collate_fn=collator,
    )
class ProtT5Cons:
    def __init__(self, checkpoint_path):
        self.model = ConsCNN().to(DEVICE)
        self.predictor = self.load_checkpoint(checkpoint_path)

    def load_checkpoint(self, checkpoint_path: Path):
        state = torch.load(checkpoint_path, map_location=DEVICE)
        self.model.load_state_dict(state["state_dict"])
        self.model = self.model.eval()
        return self.model

    def conservation_prediction(
        self, data_loader, prob_return=True, class_return=True
    ) -> ConservationOut:
        if prob_return == class_return == False:
            raise RuntimeError(
                "Method needs to either return probabilities or conservation class or both."
            )

        probability_out = dict() if prob_return else None
        class_out = dict() if class_return else None
        ret = ConservationOut(probability_out=probability_out, class_out=class_out)

        for i, (pdb_ids, X, lens) in enumerate(data_loader):
            X = X.to(DEVICE)
            with torch.no_grad():
                Yhat = self.predictor(X)
                if prob_return:
                    prob = self.model.extract_probabilities(Yhat).detach().cpu().numpy()
                if class_return:
                    cls = (
                        self.model.extract_conservation_score(Yhat)
                        .detach()
                        .cpu()
                        .numpy()
                    )
            for idx, pdb_id in enumerate(pdb_ids):
                if prob_return:
                    ret.probability_out[pdb_id] = prob[idx, :, : lens[idx]]
                if class_return:
                    ret.class_out[pdb_id] = cls[idx, : lens[idx]]
        return ret

    @staticmethod
    def write_cons_class_pred(predictions: ConservationOut, out_path):
        if predictions.class_out is None:
            raise RuntimeError("Class output not present!")
        with open(out_path, "w+") as out_f:
            for seq_id, cons_cls in predictions.class_out.items():
                out_f.write(f">{seq_id}\n{','.join([str(j) for j in cons_cls])}\n")

    @staticmethod
    def write_probabilities(predictions: ConservationOut, out_path):
        if predictions.probability_out is None:
            raise RuntimeError("Probability output not present!")
        with h5py.File(str(out_path), "w") as hf:
            for sequence_id, embedding in predictions.probability_out.items():
                # noinspection PyUnboundLocalVariable
                hf.create_dataset(sequence_id, data=embedding)
        return None
