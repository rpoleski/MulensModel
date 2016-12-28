

class LikeEvent(object):
    def __init__(self, ra=None, datasets=None):
        if isinstance(datasets, (list, tuple, LikeMulensData)) or datasets is None:
            self._set_datasets(new_value=datasets)
        else:
            raise TypeError('incorrect argument datasets of class LikeEvent()')
        if ra is not None:
            self.set_ra(ra=ra)

    @property
    def ra(self):
        return self._ra

    def set_ra(self, ra=None):
        self._ra = ra
        for dataset in self.datasets:
            dataset.set_ra_no_dependencies(ra=self._ra)

    def set_ra_no_dependencies(self, ra=None):
        self._ra = ra

    @property
    def datasets(self):
        return self._datasets

    @datasets.setter
    def datasets(self, new_value):
        self._set_datasets(new_value)

    def _set_datasets(self, new_value):
        if isinstance(new_value, LikeMulensData):
            new_value = [new_value]
        if new_value is None:
            self._datasets = None
            return
        self._datasets = new_value
        for dataset in self._datasets:
            if dataset.ra is not None:
                self.set_ra_no_dependencies(dataset.ra)
        for dataset in self._datasets:
            dataset.set_ra_no_dependencies(self.ra)

class LikeMulensData(object):
    def __init__(self, ra=None, event=None):
        if isinstance(event, LikeEvent) or event is None:
            self._set_event(event)
        else:
            raise TypeError('incorrect argument event of class LikeMulensData()')
        self.set_ra(ra=ra)

    @property
    def ra(self):
        return self._ra

    def set_ra(self, ra=None):
        self._ra = ra
        if isinstance(self.event, LikeEvent):
            self.event.set_ra_no_dependencies(ra=self._ra)

    def set_ra_no_dependencies(self, ra=None):
        self._ra = ra

    @property
    def event(self):
        return self._event

    @event.setter
    def event(self, new_value):
        self._set_event(new_value)

    def _set_event(self, new_value):
        self._event = new_value


if __name__ == '__main__':

    dat_1 = LikeMulensData(ra=10.)
    print('dat_1', dat_1.ra)
    print()

    ev_1 = LikeEvent(datasets=dat_1)
    print('ev_1', ev_1.ra)
    print('dat_1', dat_1.ra)
    print()

    dat_2 = LikeMulensData(ra=15.)
    print('dat_2', dat_2.ra)
    print()

    ev_1.datasets = [dat_1, dat_2]
    print('ev_1', ev_1.ra)
    print('dat_1', dat_1.ra)
    print('dat_2', dat_2.ra)
    print()

