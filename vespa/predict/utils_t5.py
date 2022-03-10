from transformers import T5ForConditionalGeneration, T5EncoderModel, T5Tokenizer
from transformers import logging
logging.set_verbosity_error()
from vespa.predict.config import (
    TRANSFORMER_LINK,
    DEVICE,
    EMBED, EMBEDDING_HALF_PREC, LOGODDS
)


class ProtT5:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir

    def get_model(self, model_usage: EMBED | LOGODDS):
        if model_usage == EMBED:
            model = T5EncoderModel.from_pretrained(
                TRANSFORMER_LINK, cache_dir=self.cache_dir
            )
            if EMBEDDING_HALF_PREC:
                model = model.half()
        elif model_usage == LOGODDS:
            model = T5ForConditionalGeneration.from_pretrained(
                TRANSFORMER_LINK, cache_dir=self.cache_dir
            )
        else:
            raise NotImplementedError(
                "The intended use of ProtT5 is not implemented."
            )
        model = model.eval()
        model = model.to(DEVICE)
        return model, self.get_tokenizer()

    def get_tokenizer(self):
        tokenizer = T5Tokenizer.from_pretrained(
            TRANSFORMER_LINK, do_lower_case=False, cache_dir=self.cache_dir
        )
        return tokenizer