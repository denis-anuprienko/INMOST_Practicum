import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xscale('log')
ax.set_yscale('log')

diam = [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]
errL2 = [2.22282, 0.507672, 0.12903, 0.0323996, 0.00810918, 0.00202788]
ord1 = diam
ord2 = [i ** 2 for i in diam]

ax.plot(diam, errL2, '-bo', label='$||u-u_N||_{L_2}$')
ax.plot(diam, ord1, '-r', label = '1-й порядок')
ax.plot(diam, ord2, '-g', label = '2-й порядок')
ax.legend()
plt.xlabel("Диаметр сетки")
plt.title("Для $u(x,y) = sin(10x)\cdot sin(10y) + 100x\cdot e^y$")
plt.show()
