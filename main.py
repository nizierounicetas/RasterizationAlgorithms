import time
import tkinter as tk
import matplotlib
from tkinter import messagebox

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

matplotlib.use('Agg')

class RasterizationApp:

    def __init__(self):
        self.master = tk.Tk()
        self.master.title('Rasterization algorithms')
        self.master.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.master.resizable(width=False, height=False)

        self.scale = range(0, 100)
        self.image = np.zeros((self.scale.stop, self.scale.stop))

        self.figure = plt.figure(figsize=(7, 7), dpi=100)
        plt.grid(linewidth=0.5, which='both', axis='both', color='gray', linestyle='--')
        plt.ylabel('y')
        plt.xlabel('x')
        plt.xticks(np.arange(0, self.scale.stop + 10, 10))
        plt.yticks(np.arange(0, self.scale.stop + 10, 10))
        plt.axis('equal')

        self.canvas = FigureCanvasTkAgg(self.figure)
        self.canvas.get_tk_widget().pack(side=tk.LEFT)

        self.options_frame = tk.Frame(self.master)
        self.options_frame.pack(side=tk.RIGHT, padx=10)

        self.scale_var = tk.IntVar(self.options_frame)
        self.scale_var.set(100)

        tk.Label(self.options_frame, text='Scale:').pack()
        tk.Entry(self.options_frame, textvariable=self.scale_var).pack()
        tk.Button(self.options_frame, text='Apply', command=self.change_scale).pack()

        self.algorithm_var = tk.StringVar(self.options_frame)
        self.algorithm_var.set('Step')
        algorithm_menu = tk.OptionMenu(self.options_frame, self.algorithm_var, 'Step', 'DDA', 'Bresenham',
                                       'Bresenham(circle)')
        algorithm_menu.pack()

        self.x1_var = tk.DoubleVar(self.options_frame)
        self.x1_var.set(10)
        self.y1_var = tk.DoubleVar(self.options_frame)
        self.y1_var.set(10)
        self.x2_var = tk.DoubleVar(self.options_frame)
        self.x2_var.set(40)
        self.y2_var = tk.DoubleVar(self.options_frame)
        self.y2_var.set(40)
        self.r_var = tk.DoubleVar(self.options_frame)
        self.r_var.set(50)

        tk.Label(self.options_frame, text='x1(xc):').pack()
        tk.Entry(self.options_frame, textvariable=self.x1_var).pack()
        tk.Label(self.options_frame, text='y1(yc):').pack()
        tk.Entry(self.options_frame, textvariable=self.y1_var).pack()
        tk.Label(self.options_frame, text='x2:').pack()
        tk.Entry(self.options_frame, textvariable=self.x2_var).pack()
        tk.Label(self.options_frame, text='y2:').pack()
        tk.Entry(self.options_frame, textvariable=self.y2_var).pack()
        tk.Label(self.options_frame, text='R:').pack()
        tk.Entry(self.options_frame, textvariable=self.r_var).pack()

        tk.Button(self.options_frame, text='Build', command=self.draw).pack()

        self.time_lbl = tk.Label(self.options_frame, text='Time spent:')
        self.time_lbl.pack()

        self.change_scale()


    def start(self):
        self.master.mainloop()

    def on_closing(self):
        if messagebox.askokcancel("Exit", "Are you sure you want to exit?"):
            self.master.destroy()
            plt.close()

    # Функция для пошаговой растеризации отрезка
    def rasterize_step(self):

        x1 = x2 = y1 = y2 = 0
        try:
            x1 = self.x1_var.get()
            x2 = self.x2_var.get()
            y1 = self.y1_var.get()
            y2 = self.y2_var.get()
        except:
            messagebox.showerror(message='You should enter double values', title='Error')
            return

        if not self.scale.start <= x1 <= self.scale.stop:
            messagebox.showerror(message='x1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= x2 <= self.scale.stop:
            messagebox.showerror(message='x2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y1 <= self.scale.stop:
            messagebox.showerror(message='y1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y2 <= self.scale.stop:
            messagebox.showerror(message='y2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if round(x1) == round(x2) and round(y2) == round(y1):
            self.image[round(y1), round(x1)] = 1

            plt.title('Step Algorithm')
            plt.imshow(self.image, cmap='Oranges', origin='lower')
            self.canvas.draw()
            return

        if x2 < x1:
            x2, x1 = x1, x2
            y2, y1 = y1, y2

        k = (y2 - y1) / max(x2 - x1, 0.1)
        b = y2 - k * x2
        dx = abs(x2 - x1) / (max(abs(x2 - x1), abs(y2 - y1) * 2))

        x = x1
        y = k * x + b
        self.image = np.zeros((self.scale.stop, self.scale.stop))
        begin = time.time()

        while x < x2:
            self.image[int(y), int(x)] = 1
            y = k * x + b
            x = x + dx

        end = time.time()

        plt.title('Step algorithm')
        self.time_lbl.config(text='Time spent: {}'.format(end - begin))
        plt.imshow(self.image, cmap='Oranges', origin='lower')
        self.canvas.draw()

    # Функция для растеризации отрезка по алгоритму ЦДА
    def rasterize_dda(self):

        x1 = x2 = y1 = y2 = 0
        try:
            x1 = self.x1_var.get()
            x2 = self.x2_var.get()
            y1 = self.y1_var.get()
            y2 = self.y2_var.get()
        except:
            messagebox.showerror(message='You should enter double values', title='Error')
            return

        if not self.scale.start <= x1 <= self.scale.stop:
            messagebox.showerror(message='x1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= x2 <= self.scale.stop:
            messagebox.showerror(message='x2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y1 <= self.scale.stop:
            messagebox.showerror(message='y1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y2 <= self.scale.stop:
            messagebox.showerror(message='y2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if round(x1) == round(x2) and round(y2) == round(y1):
            self.image[round(y1), round(x1)] = 1

            plt.title('Step Algorithm')
            plt.imshow(self.image, cmap='Oranges', origin='lower')
            self.canvas.draw()
            return

        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        steps = max(dx, dy)
        x_step = (x2 - x1) / steps
        y_step = (y2 - y1) / steps
        x, y = x1, y1

        self.image = np.zeros((self.scale.stop, self.scale.stop))
        begin = time.time()

        for i in range(int(steps)):
            self.image[round(y), round(x)] = 1
            x += x_step
            y += y_step

        end = time.time()

        plt.title('DDA Algorithm')
        self.time_lbl.config(text='Time spent: {}'.format(end - begin))
        plt.imshow(self.image, cmap='Oranges', origin='lower')
        self.canvas.draw()

    # Функция для растеризации отрезка по алгоритму Брезенхема
    def rasterize_bresenham(self):
        x1 = x2 = y1 = y2 = 0
        try:
            x1 = self.x1_var.get()
            x2 = self.x2_var.get()
            y1 = self.y1_var.get()
            y2 = self.y2_var.get()
        except:
            messagebox.showerror(message='You should enter double values', title='Error')
            return

        if not self.scale.start <= x1 <= self.scale.stop:
            messagebox.showerror(message='x1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= x2 <= self.scale.stop:
            messagebox.showerror(message='x2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y1 <= self.scale.stop:
            messagebox.showerror(message='y1 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= y2 <= self.scale.stop:
            messagebox.showerror(message='y2 out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        begin = time.time()

        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        sx = 1 if x1 < x2 else -1
        sy = 1 if y1 < y2 else -1
        err = dx - dy
        self.image = np.zeros((self.scale.stop, self.scale.stop))
        while True:
            if round(x1) == round(x2) and round(y1) == round(y2):
                break
            self.image[round(y1), round(x1)] = 1

            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x1 += sx
            if e2 < dx:
                err += dx
                y1 += sy

        end = time.time()

        plt.title('Bresenham algorithm')
        self.time_lbl.config(text='Time spent: {}'.format(end - begin))
        plt.imshow(self.image, cmap='Oranges', origin='lower')
        self.canvas.draw()

    # Функция для растеризации окружности по алгоритму Брезенхема
    def rasterize_circle_bresenham(self):
        xc = yc = r = 0
        try:
            xc = self.x1_var.get()
            yc = self.y1_var.get()
            r = self.r_var.get()
        except:
            messagebox.showerror(message='You should enter double values', title='Error')
            return

        if not self.scale.start <= xc <= self.scale.stop:
            messagebox.showerror(message='xc out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= yc <= self.scale.stop:
            messagebox.showerror(message='yc out of range({},{})'
                                 .format(self.scale.start, self.scale.stop), title='Error')
            return

        if not self.scale.start <= r <= (self.scale.stop / 2):
            messagebox.showerror(message='R out of range({},{})'
                                 .format(self.scale.start, (self.scale.stop / 2)), title='Error')
            return

        x, y = 0, r
        d = 3 - 2 * r

        begin = time.time()

        self.image = np.zeros((self.scale.stop, self.scale.stop))
        while x <= y:
            self.draw_circle_points(int(round(xc)), int(round(yc)), int(round(x)), int(round(y)))
            if d < 0:
                d = d + 4 * x + 6
            else:
                d = d + 4 * (x - y) + 10
                y -= 1
            x += 1

        end = time.time()

        plt.title('Bresenham algorithm(circle)')
        self.time_lbl.config(text='Time spent:{}'.format(end - begin))
        plt.imshow(self.image, cmap='Oranges', origin='lower')
        self.canvas.draw()


    # Функция для рисования точек окружности
    def draw_circle_points(self, xc, yc, x, y):
        self.image[yc + y, xc + x] = 1

        if xc - x >= 0:
            self.image[yc + y, xc - x] = 1

        if yc - y >= 0:
            self.image[yc - y, xc + x] = 1

        if yc - y >= 0 and xc - x >= 0:
            self.image[yc - y, xc - x] = 1

        self.image[yc + x, xc + y] = 1

        if xc - y >= 0:
            self.image[yc + x, xc - y] = 1

        if yc - x >= 0:
            self.image[yc - x, xc + y] = 1

        if yc - x >= 0 and xc - y >= 0:
            self.image[yc - x, xc - y] = 1

    # Функция для выбора алгоритма и настройки параметров растеризации
    def draw(self):
        # algorithm = self.algorithm_var.get()
        # if algorithm == 'Step':
        #     self.rasterize_step()
        # elif algorithm == 'DDA':
        #     self.rasterize_dda()
        # elif algorithm == 'Bresenham':
        #     self.rasterize_bresenham()
        # elif algorithm == 'Bresenham(circle)':
        #     self.rasterize_circle_bresenham()
        try:
            algorithm = self.algorithm_var.get()
            if algorithm == 'Step':
                self.rasterize_step()
            elif algorithm == 'DDA':
                self.rasterize_dda()
            elif algorithm == 'Bresenham':
                self.rasterize_bresenham()
            elif algorithm == 'Bresenham(circle)':
                self.rasterize_circle_bresenham()
        except:
            messagebox.showerror(title='Error', message='Try to increase scale')

    def change_scale(self):
        try:
            new_scale = self.scale_var.get()
            if new_scale <= 0:
                raise Exception

            self.scale = range(0, new_scale)

            plt.xticks(np.arange(0, self.scale.stop + 10, 10))
            plt.yticks(np.arange(0, self.scale.stop + 10, 10))
            self.image = np.zeros((self.scale.stop, self.scale.stop))
            plt.imshow(self.image, cmap='Oranges', origin='lower')
            self.canvas.draw()
        except:
            messagebox.showerror(message='Scale should be positive integer', title='Error')
            return

if __name__ == "__main__":
    RasterizationApp().start()


# pyinstaller --windowed -F -d bootloader main.py --name rasterization_algorithms --onefile