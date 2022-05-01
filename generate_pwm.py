import random
import math

from scipy.optimize import fsolve


def generate_single_position(ic):
    # p_1 = 0.5
    # p_a = 0.25
    # print(p_1, p_a)
    # p_t = p_1 - p_a
    # p_2 = 1 - p_1

    ic_2_p_large = {1.0:0.8105, 1.5:0.9245, 2.0:1.0}
    ic_2_p_small = {1.0:(1-0.8105)/3, 1.5:(1-0.9245)/3, 2.0:0.0}

    r = random.randint(0,3)

    result = [0.0] * 4
    for i in range(4):
        result[i] = ic_2_p_small[ic]
    result[r] = ic_2_p_large[ic]
    return result

def generate_pwm(ic, ml):
    ics = [1.0, 1.5, 2.0]
    total =  ic * ml
    pwm = []
    for col in range(ml):
        avg = total / (ml - col)
        print(avg)
        if avg == 1.0:
            new = generate_single_position(1.0)
            pwm.append(new)
            total -= 1.0
            continue
        elif avg == 2.0:
            new = generate_single_position(2.0)
            pwm.append(new)
            continue
        elif avg == 1.5:
            r = random.randint(0,2)
            new = generate_single_position(ics[r])
            pwm.append(new)
            total -= ics[r]
        elif avg > 1 and avg < 1.5:
            r = random.randint(1, 2)
            new = generate_single_position(ics[r])
            pwm.append(new)
            total -= ics[r]
            continue
        else:
            r = random.randint(0, 1)
            new = generate_single_position(ics[r])
            pwm.append(new)
            total -= ics[r]
    #print(ic)
    #print(pwm)

    return pwm

#generate_pwm(1.5, 5)
# generate_pwm(1.0, 10)
# generate_pwm(2.0, 10)
