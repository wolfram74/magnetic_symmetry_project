import numpy
from high_hepta import high_hepta_states
from low_hepta import low_hepta_states
from matplotlib import pyplot


def load_values():
    points = 10
    alpha_one_regime = []
    low_omeg_2 = []
    alpha_two_regime = []
    high_omeg_2 = []
    for i in range(points):
        alpha_one_regime.append(
            high_hepta_states[i+1][0]
            )
        low_omeg_2.append(
            high_hepta_states[i+1][2][1][0]
            )
        alpha_two_regime.append(
            low_hepta_states[-(i+1)][0]
            )
        high_omeg_2.append(
            low_hepta_states[-(i+1)][2][1][0]
            )
    return [
        numpy.array(alpha_one_regime), 
        numpy.array(low_omeg_2),
        numpy.array(alpha_two_regime),
        numpy.array(high_omeg_2)
        ]

def spot_check(omegs, alphas):
    for i in range(len(alphas)):
        print(alphas[i], omegs[i])
        if i>0:
            delta_alpha = alphas[i]-alphas[i-1]
            delta_omeg = omegs[i]-omegs[i-1] 
            print(delta_alpha, delta_omeg, delta_omeg/delta_alpha)

def slope_inter(mx_pb):
    return (mx_pb[0], -mx_pb[1]/mx_pb[0])

def relative_error(y_exp, x_exp, model):
    y_model = model[0]*x_exp+model[1]
    print(slope_inter(model))
    # print(y_model)
    # print(y_exp)
    delta = y_exp-y_model
    avg = (y_exp+y_model)#/2
    # abs_delta = abs(delta)
    abs_rel = abs(delta/avg)
    # return abs_delta/y_exp
    return abs_rel

def main():
    alpha_1, low_omeg2, alpha_2, high_omeg2 = load_values()
    # alpha_1_fit = numpy.polyfit(alpha_1, low_omeg2,  1)
    alpha_1_fit = numpy.polyfit(alpha_1, low_omeg2*alpha_1,  1)
    alpha_1_inter = slope_inter(alpha_1_fit)

    alpha_2_fit = numpy.polyfit(alpha_2, high_omeg2**2, 1)

    print(alpha_1_inter)
    rel_error1 = relative_error(low_omeg2*alpha_1, alpha_1, alpha_1_fit)
    # print(rel_error1, low_omeg2[0])
    # rel_error2 = relative_error(high_omeg2**2, alpha_2, alpha_2_fit)
    # print(rel_error2)

    # print(slope_inter(alpha_1_fit))

    # plotting(low_omeg2, alpha_1, alpha_1_fit)
    # plotting(low_omeg2*alpha_1, alpha_1, alpha_1_fit)

    # plotting(high_omeg2, alpha_2, alpha_2_fit, sqrt=True)
    # plotting(high_omeg2**2, alpha_2, alpha_2_fit, power =4)

def plotting(y_exp, x_exp, model, sqrt=False, power=2):

    figure, axes = pyplot.subplots()
    fontsize=20
    axes.set_xlabel('$\\alpha$', fontsize=fontsize)
    axes.set_ylabel('$(\\omega_2/\\Omega)^%d$' % power, fontsize=fontsize)
    axes.tick_params(axis='both', labelsize=14)


    x_model = x_exp
    if power==2:
        x_model = list(x_model)
        x_model.insert(0, 1.15)
        x_model=numpy.array(x_model)

    y_model = model[0]*x_model+model[1]
    if sqrt:
        y_model = numpy.sqrt(y_model)

    axes.plot(x_exp, y_exp, "x")
    axes.plot(x_model, y_model)
    axes.axhline(y=0,color='k',linestyle='--')
    if power==2:
        yticks = (0, .05, .1, .15, .2)
        axes.set_yticks(yticks)

    pyplot.tight_layout()
    pyplot.savefig("comparison_power_%d_mk4.pdf" % power)

if __name__ == '__main__':
    main()

'''
# strogatz problem 3.4.16

fri donburri
sat salmon
sun pizza
mon quiche
tue grill

'''