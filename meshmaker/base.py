class Base:

    def __init__(self, **kws):
        for k, v in kws.items():
            if not hasattr(self, k):
                setattr(self, k, v)
