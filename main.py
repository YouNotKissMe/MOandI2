from scipy.optimize import golden
from tabulate import tabulate
from math import e, pi, sqrt,fabs
from sympy import *
# Найти локальный минимум функции: f(X) = x1^2 + (x2 - 3)^2 x0 = (0,5 ; 4) eps = 0,01
# a) метод покоординатного спуска(метод конфигруаций)
# б) метод наискорейшего спуска
# в) метод сопреженных гардиентов


class Gardient:
    def __init__(self, start: list, eps):
        self.x,self.y = Symbol('x'),Symbol('y')
        self.diffx, self.diffy = diff(2* self.x ** 2 + self.x*self.y+ self.y ** 2, self.x), \
                                 diff(2* self.x ** 2 + self.x*self.y+ self.y ** 2,self.y)  # здесь находим производную от функции
        self.__x1,self.__x2 = start[0], start[1]
        self.__startX1, self.__startX2 = antigardient(self.diffx,  self.__x1, self.__x2), \
                                         antigardient(self.diffy,  self.__x1, self.__x2)
        self.__eps = eps
        self.__k = 0
        self.b1, self.b2 = None , None
        self.table = [['N', 'x(n)', 'S(n)'], [self.__k, f'({self.__x1};{self.__x2})', f'({self.__startX1};{self.__startX2})']]

    def algoritm(self):   #метод наискорейшего спуска
        if sqrt(self.__startX2 ** 2 + self.__startX1 ** 2) < self.__eps:
            print(sqrt(self.__startX2 ** 2 + self.__startX1 ** 2))
            print(tabulate(self.table, tablefmt='pipe', stralign='center',
                           headers='firstrow'))
        else:
            self.__k += 1

            def func(x):
               #лямда функция, которая приходит на замену 2 переменным
                return 2 * (self.__x1 + self.__startX1 * x)**2 + \
                       (self.__x1 + self.__startX1 * x)*(self.__x2 + self.__startX2 * x) + \
                       (self.__x2 + self.__startX2 * x)** 2
            minimize = golden(func)
            self.__x1, self.__x2 = self.__x1 + self.__startX1 *\
                                   minimize,\
                                   self.__x2 + self.__startX2 *\
                                   minimize
            self.__startX1, self.__startX2 = antigardient(self.diffx,  self.__x1, self.__x2), \
                                            antigardient(self.diffy,  self.__x1, self.__x2)
            self.table.append([self.__k, f'({self.__x1};{self.__x2})', f'({self.__startX1};{self.__startX2})'])
            Gardient.algoritm(self)

    def rivz(self):
        if sqrt(self.__startX2 ** 2 + self.__startX1 ** 2) < self.__eps:
            print(sqrt(self.__startX2 ** 2 + self.__startX1 ** 2))
            print(tabulate(self.table, tablefmt='pipe', stralign='center',
                           headers='firstrow'))
        else:
            self.__k += 1

            def func(x):
                # лямда функция, которая приходит на замену 2 переменным
                return 2 * (self.__x1 + self.__startX1 * x) ** 2 + \
                       (self.__x1 + self.__startX1 * x) * (self.__x2 + self.__startX2 * x) + \
                       (self.__x2 + self.__startX2 * x) ** 2

            minimize = golden(func)
            self.__x1, self.__x2 = self.__x1 + self.__startX1 * \
                                   minimize, \
                                   self.__x2 + self.__startX2 * \
                                   minimize
            self.b1, self.b2 = (antigardient(self.diffx, self.__x1, self.__x2)**2 / self.__startX1**2) ,\
                               (antigardient(self.diffy, self.__x1, self.__x2**2)  /self.__startX2**2)
            self.__startX1, self.__startX2 = antigardient(self.diffx, self.__x1, self.__x2) + self.b1 * self.__startX1, \
                                             antigardient(self.diffy, self.__x1, self.__x2)+ self.b2 * self.__startX2
            self.table.append([self.__k, f'({self.__x1};{self.__x2})', f'({self.__startX1};{self.__startX2})'])
            Gardient.rivz(self)


def antigardient(x,x2,y):
    x = str(x).split()
    m = ''
    for i in x:
        m += i
    m1 = ''
    for i in m:
        if i == 'x':
            m1 += str(x2)
        elif i == 'y':
            m1 += str(y)
        else:
            m1 += i
    return -1 * (eval(m1))


Gardient([0.5, 1], 0.1).algoritm()
Gardient([0.5, 1], 0.1).rivz()