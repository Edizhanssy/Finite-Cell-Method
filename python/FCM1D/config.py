import numpy as np
from dataclasses import dataclass, field

@dataclass(frozen=True)
class Config:
    L: float = 3.0
    n_elements: int = 2
    E: float = 1.0
    A: float = 1.0
    alpha: float = 1e-8
    p: int = 4
    max_depth: int = 4
    penalty: float = 1e5

    # (x_start, x_end, material_id) - 0: physical, 1: fictitious
    domains: tuple = ((0.0, 1.0, 0), (1.0, 7/3, 1), (7/3, 3.0, 0))

    load_span: tuple = (0.0, 1.0)
    load_amp: float = 1/20
    load_freq: float = 4 * np.pi


    @property
    def n_gauss(self):
        return self.p + 1
    @property
    def n_modes(self):
        return self.p + 1

    def domain_index(self, x):
        for x0, x1, mat in self.domains:
            if x0 <= x <= x1:
                return mat
        return -1

    def material_factor(self, material_id):
        if material_id == 0:
            return 1.0
        if material_id == 1:
            return self.alpha
        raise ValueError(f'unvalid material id: {material_id}')

    def body_load(self, x):
        a, b = self.load_span
        return self.load_amp * np.sin(self.load_freq * x) if a <= x <= b else 0.0

    disp_load: float = -1.0

    def to_local(self, x):
        return (2 * x - self.L) / self.L

    @property
    def fictitious_span(self):
        for x0, x1, mat in self.domains:
            if mat == 1:
                return (x0, x1)
        raise ValueError('fictitious domain is not defined')