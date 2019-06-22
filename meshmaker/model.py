from collections import defaultdict


class model:

    counter = 0

    def __new__(cls):
        name = f'model_{cls.counter}'
        cls.counter += 1
        inst = super().__new__(cls)
        inst.__init__()
        inst.name = name
        print('name', inst.name)
        return inst

    def __init__(self):
        self.surfaces = defaultdict(list)

    def add(self, surface, material):
        self.surfaces[material].append(surface)
        return self


