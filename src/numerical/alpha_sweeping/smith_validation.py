import system
import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()


def get_vector_mu():
    output = []
    magnets.load_state(.01)
    while magnets.alpha <2.0:
        values = [magnets.alpha]
        values.append([(cos(el), sin(el)) for el in magnets.state[:7]])
        output.append(values)
        magnets.load_state(magnets.alpha+0.011)
    for values in output[:5]:
        print(values)
    return output

def validate_symmetries():
    vectorized_mus = get_vector_mu()
    template = "%.4f, "+"%.4f, "*6
    print('alpha, m1x-m4x, m2x-m3x, m5x-m6x, m1y+m4y, m2y+m3y, m5y+m6y')
    for vals in vectorized_mus:
        comparisons = [vals[0],
        vals[1][1][0]-vals[1][4][0],vals[1][2][0]-vals[1][3][0],vals[1][5][0]-vals[1][6][0],
        vals[1][1][1]+vals[1][4][1],vals[1][2][1]+vals[1][3][1],vals[1][5][1]+vals[1][6][1]]
        print(template%tuple(comparisons))

if __name__ == '__main__':
    # get_vector_mu()
    validate_symmetries()

