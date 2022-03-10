import numpy
import torch
import h5py
from tqdm import tqdm
from pathlib import Path

from vespa.predict.config import (
    DEVICE, CACHE_DIR, VERBOSE,
    EMBED, EMB_MAX_SEQ_LEN, EMB_MAX_RESIDUES, EMB_MAX_BATCH, EMB_STORE_FREQ
)
from vespa.predict.utils import parse_fasta_input
from vespa.predict.utils_t5 import ProtT5


class T5_Embed:
    def __init__(self, cache_dir):
        self.prott5 = ProtT5(cache_dir)
        self.saving_pattern = 'w'

    def embed_from_fasta(self, fasta_path, output_path):
        self.saving_pattern = 'w'
        if VERBOSE:
            print('Load model: ProtT5')
        self.model, self.tokenizer = self.prott5.get_model(EMBED)
        if VERBOSE:
            print('Compute embeddings!')
        self.get_embeddings(fasta_path, output_path)

    def embedding_init(self, fasta_path):
        seq_dict = parse_fasta_input(fasta_path)
        seq_dict = sorted(seq_dict.items(), key=lambda kv: len(seq_dict[kv[0]]), reverse=True)
        return seq_dict

    def process_batch(self, batch, emb_dict):
        pdb_ids, seqs, seq_lens = zip(*batch)

        token_encoding = self.tokenizer(seqs, add_special_tokens=True, padding='longest', return_tensors="pt")
        input_ids = token_encoding['input_ids'].to(DEVICE)
        attention_mask = token_encoding['attention_mask'].to(DEVICE)

        try:
            # batch-size x seq_len x embedding_dim
            with torch.no_grad():
                embedding_repr = self.model(input_ids, attention_mask=attention_mask)
        except RuntimeError:
            print("RuntimeError for {} (L={})".format(pdb_ids, seq_lens))
            return emb_dict

        new_emb_dict = dict()
        for batch_idx, identifier in enumerate(pdb_ids):
            s_len = seq_lens[batch_idx]
            emb = embedding_repr.last_hidden_state[batch_idx, :s_len]
            new_emb_dict[identifier] = emb.detach().cpu().numpy().squeeze()

        if new_emb_dict:
            emb_dict.update(new_emb_dict)
        return emb_dict

    def save_embeddings(self, output_path, emb_dict):
        Path(str(output_path.absolute())).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(str(output_path.absolute()), self.saving_pattern) as hf:
            for sequence_id, embedding in emb_dict.items():
                hf.create_dataset(sequence_id, data=embedding)
        self.saving_pattern = 'a'

    def get_embeddings(self, fasta_path, output_path):
        seq_dict = self.embedding_init(fasta_path)

        emb_dict = dict()
        batch, n_res_batch = [], 0

        for seq_idx, (pdb_id, seq) in tqdm(enumerate(seq_dict, 1), total=len(seq_dict)):
            seq_len = len(seq)
            seq = ' '.join(list(seq))

            if seq_len >= EMB_MAX_SEQ_LEN:
                emb_dict = self.process_batch([(pdb_id, seq, seq_len)], emb_dict)
            else:
                if len(batch) >= EMB_MAX_BATCH or n_res_batch >= EMB_MAX_RESIDUES:
                    emb_dict = self.process_batch(batch, emb_dict)
                    batch = []
                    n_res_batch = 0

                batch.append((pdb_id, seq, seq_len))
                n_res_batch += seq_len

            if len(emb_dict) > EMB_STORE_FREQ:
                self.save_embeddings(output_path, emb_dict)
                emb_dict = dict()

        if batch:
            emb_dict = self.process_batch(batch, emb_dict)

        if emb_dict:
            self.save_embeddings(output_path, emb_dict)
