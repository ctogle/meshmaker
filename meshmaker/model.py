from collections import defaultdict


class model:

    counter = 0

    def __init__(self):
        self.surfaces = defaultdict(list)
        self.name = 'model_{}'.format(model.counter)
        model.counter += 1

    def add(self, surface, material):
        self.surfaces[material].append(surface)
        return self


