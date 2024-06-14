import math

# CONFIG
STEP = 0.001


def calc_theta(p1, p2, p3):
    hypotenuse1 = math.dist(p2, p1)
    height1 = p2[1] - p1[1]  # y coord
    angle1 = math.pi - math.asin(height1 / hypotenuse1)
    
    hypotenuse2 = math.dist(p3, p2)
    height2 = p3[1] - p2[1]  # y coord
    angle2 = math.pi - math.asin(height2 / hypotenuse2)
    
    return angle1 - angle2


# Turning angle
kappa_a = lambda theta: theta

# Length variation
kappa_b = lambda theta: 2 * math.sin(theta / 2)

# Steiner formula
kappa_c = lambda theta: 2 * math.tan(theta / 2)

# Osculating circle
kappa_d = lambda theta, w: 2 * math.sin(theta) / w


def calc_kappa(step):
    samples = [0, step, 2*step]
    curve_section = [(math.cos(t), math.sin(t)) for t in samples]

    theta = -1 * calc_theta(*curve_section)
    w = math.dist(curve_section[-1], curve_section[0])
    
    return [kappa_a(theta), kappa_b(theta), kappa_c(theta), kappa_d(theta, w)]


def main():
    print(calc_kappa(STEP))


if __name__ == '__main__':
    main()
