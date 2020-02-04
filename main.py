from kivy.app import App
from kivy.uix.widget import Widget
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.popup import Popup
from kivy.uix.textinput import TextInput
from kivy.properties import ListProperty
import numpy as np
from random import shuffle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from math import *
import threading

class GenWidget(Widget):
    result = ListProperty()

class ind:
    indCount = 0
    rank = 0

    def __init__(self, genes):
        self.genes = genes
        self.rank = 0
        self.fitness = 0


class GenApp(App):
    def build(self):
        return GenWidget()

    def SUS(self, Population):
        F = sum(ind.rank for ind in Population)
        # N := liczba osobników do zachowania
        N = len(Population) / 2 + 1
        # P: = distance between the pointers
        P = F / N
        Start = np.random.rand() * P
        pointers = []
        for i in range(int(N - 1)):
            pointers.append(Start + i * P)
        return (self.RWS(Population, pointers))

    def RWS(self, Population, Points):
        Keep = []
        for P in Points:
            i = 0
            while sum(ind.rank for ind in Population[:i]) < P:
                i = i + 1
            Keep.append(Population[i])
        return (Keep)

    def delta(self, t, y, T, r):
        return (y * (1 - r ** ((1 - t / T) * 2)))

    def evolve(self, cel, genrange, individuals, end_n, eps_gen, eps_fen, maxstate, plot_def):
        indies = []
        genrange = list(map(float, genrange.split(", ")))
        splitind = cel.index(":")
        fargs = cel[0:splitind]
        fargs2 = fargs
        fargs = list(map(str, fargs.split(", ")))
        individuals = int(individuals)
        end_n = int(end_n)
        eps_gen = float(eps_gen)

        for i in range(individuals):
            chrom = []
            for j in range(0, int(len(genrange) / 2)):
                chrom.append((genrange[2 * j + 1] - genrange[2 * j]) * np.random.rand() - genrange[2 * j])
            indies.append(ind(chrom))

        eps_fen = float(eps_fen)

        for i in range(0, len(indies)):
            for j in range(len(fargs)):
                # setattr is used, because we don't know how many function's attributes will be
                setattr(indies[i], fargs[j], indies[i].genes[j])
            code = '(' + 'lambda ' + cel + ')(' + str((indies[i].genes)).strip('[]') + ')'
            code = code.strip('[]')
            indies[i].fitness = eval(code)

        for i in range(end_n):
            if maxstate == "down":
                indies.sort(key=lambda ind: ind.fitness, reverse=True)
            else:
                indies.sort(key=lambda ind: ind.fitness, reverse=False)
            for j in range(len(indies)):
                indies[j].rank = len(indies) - j
            parents = self.SUS(indies)

            a = 0.25
            shuffle(parents)
            indies = []
            indies = parents[:]
            chrom1 = []
            chrom2 = []
            x = int(len(parents) / 2)
            for l in range(0, x):
                l2 = l + x
                k = np.random.randint(0, len(fargs))
                chrom1 = parents[l].genes[:k]
                chrom2 = parents[int(l + len(parents) / 2)].genes[:k]
                for m in range(k, len(parents[l].genes)):
                    chrom1.append(a * parents[l].genes[m] + (1 - a) * parents[l2].genes[m])
                    chrom2.append(a * parents[l].genes[m] + (1 - a) * parents[l2].genes[m])

                # mutacja
                for t in range(len(chrom1)):
                    r = np.random.rand()
                    if np.random.rand() < 0.5:
                        chrom1[t] = chrom1[t] + self.delta(i, genrange[2 * t + 1] - chrom1[t], end_n, r)
                    else:
                        chrom1[t] = chrom1[t] - self.delta(i, chrom1[t] - genrange[2 * t], end_n, r)
                    r = np.random.rand()
                    if np.random.rand() < 0.5:
                        chrom2[t] = chrom2[t] + self.delta(i, genrange[2 * t + 1] - chrom2[t], end_n, r)
                    else:
                        chrom2[t] = chrom2[t] - self.delta(i, chrom2[t] - genrange[2 * t], end_n, r)

                indies.append(ind(chrom1))
                indies.append(ind(chrom2))
            for n in range(0, len(indies)):
                for j in range(len(fargs)):
                    setattr(indies[n], fargs[j], indies[n].genes[j])
                code = '(' + 'lambda ' + cel + ')(' + str((indies[n].genes)).strip('[]') + ')'
                code = code.strip('[]')
                setattr(indies[n], 'fitness', eval(code))

            meanfen = sum(ind.fitness for ind in indies) / len(indies)
            meangen = []
            for g in range(len(fargs)):
                meangen.append(sum(ind.genes[g] for ind in indies) / len(indies))
            diff = list(np.array(meangen)-np.array(indies[0].genes))
            diff = map(abs, diff)
            cause = "else"
            if (all(i <= eps_gen for i in diff)):
                cause = "genotype similarity"
                break
            if (abs(indies[0].fitness - meanfen) <= eps_fen):
                cause = "fenotype similarity"
                break
            if (cause != "genotype similarity" and cause != "fenotype similarity"):
                cause = "iteration limit reached"

        if maxstate == "down":
            indies.sort(key=lambda ind: ind.fitness, reverse=True)
        else:
            indies.sort(key=lambda ind: ind.fitness, reverse=False)
        result = ("[" + str(fargs2) + "]" + "=" + str(indies[0].genes) + "\n" + "fitness = " + str(indies[0].fitness) + "   end: " + cause)

        if (plot_def == "down"):
            if (len(fargs) == 1):
                x = np.linspace(genrange[0], genrange[1], 100)
                y = np.zeros(len(x))
                for n in range(0, len(x)):
                    y[n] = eval('(' + 'lambda ' + cel + ')(' + str(x[n]) + ')')
                plt.plot(x, y)
                plt.scatter(indies[0].genes[0], indies[0].fitness, edgecolors="r")
                plt.xlabel("x")
                plt.ylabel("y")
                plt.savefig((time.strftime("%d_%m_%Y")+time.strftime("_%I_%M_%S")))
                plt.close()
            if (len(fargs) == 2):
                x = np.linspace(genrange[0], genrange[1], 50)
                y = np.linspace(genrange[2], genrange[3], 50)
                z = np.zeros((len(x), len(y)))
                for n in range(len(x)):
                    for m in range(0, len(y)):
                        z[n][m] = eval('(' + 'lambda ' + cel + ')(' + str(x[n]) +  "," + str(y[m]) + ')')
                X, Y = np.meshgrid(x, y)
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_surface(X, Y, z, color='b')
                ax.scatter([indies[0].genes[1]],[indies[0].genes[0]], [indies[0].fitness], color="r")
                plt.xlabel("y")
                plt.ylabel("x")
                plt.savefig((time.strftime("%d_%m_%Y")+time.strftime("_%I_%M_%S")))
                plt.close()
        return (result)

if __name__ == '__main__':
    GenApp().run()