from inputprograms import NASACEA

class CEAInterface:
    @staticmethod
    def compute(Pc, OF, epsilon):
        return NASACEA.CEA(Pc, OF, epsilon)
