from operator import xor
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from mpl_toolkits.mplot3d import Axes3D


def get_filename():
    path = './data'
    file_list = []
    for root, dirs, files in os.walk(path):
        return files

def plot_conv_rate(content):
    data = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    for t in content[3]:
        data.append(float(t))
    if(len(data) != 4):
        return -1
    x = [5, 6, 7, 8]
    data = np.log2(np.array(data))
    plt.clf()
    plt.plot(x, data)
    plt.title(title)
    plt.xticks(x)
    plt.xlabel("$log_2(N)$")
    plt.ylabel("$log_2(||err||_\infty)$")
    plt.savefig(save_path)
    return 0

def plot_1d_error(content):
    data = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    for t in content[3]:
        data.append(float(t))
    x = []
    s = 1/float(len(data)-1)
    for i in range(len(data)):
        x.append(i*s)
    # plt.clf()
    plt.clf()
    plt.plot(x, data)
    plt.title(title)
    plt.xlabel('$x$ axis')
    plt.ylabel('absolute errors')
    plt.savefig(save_path)
    return 0

def plot_2d_error(content):
    data = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    for t in content[3]:
        data.append(float(t))
    data = np.array(data)
    N = round(np.sqrt(len(data)))
    if(N*N != len(data)):
        return -1
    data = data.reshape(N, N)
    # print(data)
    x = np.arange(0, 1+1e-8, 1/float(N-1))
    X, Y = np.meshgrid(x, x)
    fig = plt.figure()
    ax = Axes3D(fig,auto_add_to_figure=False, title = title)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, data, rstride = 1, cstride = 1, cmap = plt.get_cmap('rainbow'))
    ax.set_xlabel('$x$ axis')
    ax.set_ylabel('$y$ axis')
    ax.set_zlabel('absolute errors')
    ax.set_title(title, y = 1)
    plt.savefig(save_path)
    return 0

def plot_VC(content):
    data = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    for t in content[3]:
        data.append(float(t))
    # x = []
    # s = 1/float(len(data)-1)
    # for i in range(len(data)):
    #     x.append(i*s)
    data = np.array(data)
    plt.clf()
    plt.plot(range(len(data)), data)
    plt.title(title)
    plt.xlabel('the number of iterations')
    plt.ylabel('maximum norm of residual ($||r||_{\infty}$)')
    plt.savefig(save_path)
    plt.clf()
    plt.plot(range(len(data)), np.log10(data))
    plt.title(title)
    plt.xlabel('the number of iterations')
    plt.ylabel('$log_{10}(||r||_{\infty})$')
    plt.savefig(save_path+'log')
    # print(data)
    return 0

def plot_iter(content):
    data = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    for t in content[3]:
        data.append(float(t))
    # x = np.arange(8, 8+len(data)*0.5, 0.5)
    x = range(8, 8+len(data), 1)
    plt.clf()
    plt.plot(x, data)
    plt.title(title)
    plt.xticks(x)
    plt.xlabel('$-log_{10}(\\varepsilon)$')
    plt.ylabel('the number of iterations')
    plt.savefig(save_path)
    # print(data)
    return 0

def plot_time(content):
    VC = []
    FMG = []
    LU = []
    title = content[0][0]
    save_path = './fig/'+content[1][0]
    # 32 = 2^5
    N = int(len(content[3])/3)
    for i in range(N):
        VC.append(float(content[3][i*3]))
        FMG.append(float(content[3][i*3+1]))
        LU.append(float(content[3][i*3+2]))
    plt.clf()
    plt.plot(range(5, 5+N, 1), VC, label = 'V-cycle')
    plt.plot(range(5, 5+N, 1), FMG, label = 'FMG')
    plt.plot(range(5, 5+N, 1), LU, label = 'LU')
    plt.xlabel('$log_2(N)$')
    plt.ylabel('time (second)')
    plt.xticks(range(5,5+N,1))
    plt.title(title)
    plt.legend()
    plt.savefig(save_path)
    return 0

def plot(filename):
    path = './data/'+filename
    csv_reader = csv.reader(open(path))
    content = []
    for row in csv_reader:
        content.append(row)
    type = content[2][0]
    dim = int(content[2][1])
    if(type == "err"):
        if(dim == 1):
            plot_1d_error(content)
        else:
            plot_2d_error(content)
    elif(type == 'VC') :
        plot_VC(content)
    elif(type == 'convergence'):
        plot_conv_rate(content)
    elif(type == 'iter'):
        plot_iter(content)
    elif(type == 'time'):
        plot_time(content)
    return 0

def plot_all():
    file_list = get_filename()
    for filename in file_list:
        try:
            plot(filename)
        except:
            continue
    return 0

# plot('test1main.csv')
plot_all()

# plot('2d_dirichletCPUtime.csv')

# plot("1d_mixediter.csv")

# plot('test2_dirichletabserr_N32.csv')
# plot('test2_dirichletabserr_N64.csv')
# plot('test2_dirichletabserr_N128.csv')
# plot('test2_dirichletabserr_N256.csv')

# plot('test2_irabserr_N32.csv')
# plot('test2_irabserr_N64.csv')
# plot('test2_irabserr_N128.csv')
# plot('test2_irabserr_N256.csv')

'''
(title) xxxx
(filename) xxx.png
(type) err/VC/convergence, 1(dim)
data ...
'''