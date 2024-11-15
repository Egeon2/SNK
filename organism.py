from random import randint
from tqdm import tqdm
import matplotlib.pyplot as plt
import random
import numpy as np
from matplotlib.animation import FuncAnimation

nucleotides_dna = np.array(['A', 'T', 'G', 'C'])

class Cell:

    def __init__(self, dna, comp_factor, antibiotics, coords_x, coords_y, conc_antibiotics):
        """
        Инициализация клетки.
        :param dna: Строка с последовательностью ДНК клетки.
        :param comp_factor: Количество компетентного фактора.
        :param antibiotics: Количество молекул антибиотиков.
        :param coords_x: Координата X клетки.
        :param coords_y: Координата Y клетки.
        :param conc_antibiotics: Концентрация антибиотиков.
        """
        self.dna = dna
        self.comp_factor = comp_factor
        self.antibiotics = antibiotics
        self.coords_x = coords_x
        self.coords_y = coords_y
        self.conc_antibiotics = conc_antibiotics
        self.alive = True

    # Группа 1: Методы, связанные с состоянием клетки (жива ли она и её смерть)
    def check_alive(self, resist_gen):
        """
        Проверка, жива ли клетка, в зависимости от концентрации антибиотика и наличия гена резистентности.
        """
        if self.conc_antibiotics < 10 and resist_gen in self.dna:
            print("Cell is alive")
            self.alive = True
        elif resist_gen in self.dna and self.conc_antibiotics > 10:
            print("Cell is alive")
            self.alive = True
        else:
            print("Cell is dead")
            self.alive = False
            self.fragment_dna()  # Фрагментация ДНК при смерти
            self.drop_dna_on_death()  # Отлет фрагментов ДНК при смерти
        return self.alive

    def fragment_dna(self):
        """
        Делит ДНК клетки на случайные фрагменты, если клетка погибла.
        """
        dna_length = len(self.dna)
        fragments = []
        i = 0
        while i < dna_length:
            frag_length = np.random.randint(1, 5)  # Длина фрагмента от 1 до 5 нуклеотидов
            fragments.append(self.dna[i:i + frag_length])
            i += frag_length
        print("ДНК фрагментирована:", fragments)
        return fragments

    # Группа 2: Методы, связанные с дубликацией и трансформацией ДНК
    def duplication(self, num_of_muts):
        temp_res = list(self.dna)
        position = np.random.randint(0, len(temp_res)-1)
        temp_res[position] = np.random.choice(nucleotides_dna, num_of_muts)
        dupl_dna = temp_res
        return temp_res, dupl_dna

    def transform(self, donor_dna, donor_x, donor_y, max_distance=10):
        """
        Метод трансформации, который проверяет наличие гомологичных участков и выполняет рекомбинацию.
        :param donor_dna: Последовательность ДНК донора.
        :param donor_x: Координата X донора.
        :param donor_y: Координата Y донора.
        :param max_distance: Максимальное расстояние для трансформации.
        :return: None
        """
        distance = self.distance_to(donor_x, donor_y)
        print(f"Расстояние до донора: {distance}")
        if distance <= max_distance:
            if self.comp_factor > 0:
                is_homologous = self.check_homology(donor_dna)
                if is_homologous:
                    print("Гомологичная ДНК обнаружена. Встраивание...")
                    self.recombine(donor_dna)
                else:
                    print("Негомологичная ДНК. Выбрасывание...")
                    self.drop_dna(self.coords_x, self.coords_y)  # Выбрасываем негомологичную ДНК
            else:
                print("Недостаточно компетентного фактора для трансформации.")
        else:
            print("Слишком большое расстояние для трансформации.")

    def check_homology(self, donor_dna):
        """
        Метод проверки гомологии между ДНК реципиента и донора.
        :param donor_dna: Последовательность ДНК донора.
        :return: True, если ДНК гомологичная, иначе False.
        """
        return random.random() > 0.5  # Замените на реальную проверку гомологии

    def recombine(self, donor_dna):
        """
        Метод рекомбинации для встраивания гомологичной ДНК.
        :param donor_dna: Гомологичная последовательность ДНК донора.
        :return: None
        """
        self.dna = self.dna.replace("гомологичный_участок", donor_dna)
        print("Рекомбинация завершена. ДНК обновлена:", self.dna)

    # Группа 3: Методы, связанные с координатами и расстояниями
    def coordinates(self):
        return self.coords_x, self.coords_y

    def distance_to(self, other_x, other_y):
        """
        Вычисляет евклидово расстояние до других координат.
        :param other_x: Координата X другой клетки.
        :param other_y: Координата Y другой клетки.
        :return: Евклидово расстояние.
        """
        return np.sqrt((self.coords_x - other_x) ** 2 + (self.coords_y - other_y) ** 2)

    def drop_dna(self, drop_coords_x, drop_coords_y):
        """
        Выбрасывает ДНК клетки на случайные координаты.
        """
        drop_coords_x = np.random.randint(1, 3)
        drop_coords_y = np.random.randint(1, 3)
        print(f"ДНК выброшена на новые координаты: ({drop_coords_x}, {drop_coords_y})")
        return drop_coords_x, drop_coords_y

    def drop_dna_on_death(self):
        """
        При смерти клетки фрагменты ДНК выбрасываются на случайные координаты.
        """
        print("ДНК фрагментирована и выброшена...")
        fragments = self.fragment_dna()  # Получаем фрагменты ДНК
        for frag in fragments:
            drop_coords_x, drop_coords_y = self.drop_dna(self.coords_x, self.coords_y)
            print(f"Фрагмент '{frag}' выброшен на координаты: ({drop_coords_x}, {drop_coords_y})")

    def diffuse_antibiotics(self, grid, diffusion_rate=0.1):
        """
        Моделирование диффузии антибиотиков на сетке.
        :param grid: Двумерный массив (numpy.ndarray), представляющий концентрации антибиотиков.
        :param diffusion_rate: Коэффициент диффузии (доля перемещаемого вещества).
        :return: Обновленная сетка.
        """
        new_grid = grid.copy()
        grid_size = grid.shape[0]

        for x in range(grid_size):
            for y in range(grid_size):
                if grid[x, y] > 0:  # Если в ячейке есть антибиотик
                    diffusion_amount = grid[x, y] * diffusion_rate

                    # Распределяем часть концентрации на соседей (4 соседние ячейки)
                    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        nx, ny = x + dx, y + dy
                        if 0 <= nx < grid_size and 0 <= ny < grid_size:
                            new_grid[nx, ny] += diffusion_amount / 4

                    # Уменьшаем текущую концентрацию
                    new_grid[x, y] -= diffusion_amount

        return new_grid


# Размер сетки
grid_size = 100  # Устанавливаем размер сетки 10x10
grid = np.zeros((grid_size, grid_size))  # Инициализация пустой сетки (0 - пустое место)

# Функция для генерации случайных координат
def generate_random_coords(grid_size, count):
    return [(random.randint(0, grid_size - 1), random.randint(0, grid_size - 1)) for _ in range(count)]

# Функция для добавления клеток
def add_cells(grid, cell_count, is_recipient=False):
    cells = []
    for _ in range(cell_count):
        coords = generate_random_coords(grid_size, 1)[0]
        dna = "реципиент_ДНК" if is_recipient else f"донорная_ДНК{random.randint(1, 10)}"
        comp_factor = random.randint(0, 5) if is_recipient else 0
        antibiotics = random.randint(0, 3)
        conc_antibiotics = random.randint(0, 20) if is_recipient else 0
        cell = {"dna": dna, "comp_factor": comp_factor, "antibiotics": antibiotics, "x": coords[0], "y": coords[1], "conc_antibiotics": conc_antibiotics}
        cells.append(cell)
        grid[coords[0], coords[1]] = 2 if is_recipient else 3  # 2 - реципиент, 3 - донор
    return cells

# Функция для добавления антибиотиков
def add_antibiotics(grid, antibiotic_count):
    antibiotic_coords = generate_random_coords(grid_size, antibiotic_count)
    for x, y in antibiotic_coords:
        grid[x][y] = 1  # 1 - антибиотик
    return antibiotic_coords

# Инициализация объектов на сетке
recipient_cells = add_cells(grid, 500, is_recipient=True)
donor_cells = add_cells(grid, 500)
antibiotic_coords = add_antibiotics(grid, 500)

# Настройка визуализации
fig, ax = plt.subplots(figsize=(8, 8))
cmap = plt.get_cmap("viridis")
ax.set_title("Карта клеток и антибиотиков")
im = ax.imshow(grid, cmap=cmap, interpolation='nearest')

def update(frame):
    global grid
    new_grid = grid.copy()  # Копируем текущую сетку для обновлений

    # Один шаг диффузии для антибиотиков 1, 2 и 3
    for x in range(grid_size):
        for y in range(grid_size):
            if grid[x, y] == 1:  # Если в клетке есть антибиотик 1
                possible_moves = []  # Список возможных направлений для движения антибиотика

                # Проверка соседей (вверх, вниз, влево, вправо)
                if x - 1 >= 0 and grid[x - 1, y] == 0:
                    possible_moves.append((-1, 0))  # Вверх
                if x + 1 < grid_size and grid[x + 1, y] == 0:
                    possible_moves.append((1, 0))  # Вниз
                if y - 1 >= 0 and grid[x, y - 1] == 0:
                    possible_moves.append((0, -1))  # Влево
                if y + 1 < grid_size and grid[x, y + 1] == 0:
                    possible_moves.append((0, 1))  # Вправо

                # Если есть возможные движения, случайным образом выбрать одно
                if possible_moves:
                    dx, dy = random.choice(possible_moves)
                    nx, ny = x + dx, y + dy

                    # Дополнительная проверка на занятость клетки
                    if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                        new_grid[nx, ny] = 1  # Перемещаем антибиотик в новую клетку
                        new_grid[x, y] = 0  # Обнуляем старую клетку

            elif grid[x, y] == 2:  # Если в клетке есть антибиотик 2
                possible_moves = []  # Список возможных направлений для движения антибиотика

                # Проверка соседей (вверх, вниз, влево, вправо)
                if x - 1 >= 0 and grid[x - 1, y] == 0:
                    possible_moves.append((-1, 0))  # Вверх
                if x + 1 < grid_size and grid[x + 1, y] == 0:
                    possible_moves.append((1, 0))  # Вниз
                if y - 1 >= 0 and grid[x, y - 1] == 0:
                    possible_moves.append((0, -1))  # Влево
                if y + 1 < grid_size and grid[x, y + 1] == 0:
                    possible_moves.append((0, 1))  # Вправо

                # Если есть возможные движения, случайным образом выбрать одно
                if possible_moves:
                    dx, dy = random.choice(possible_moves)
                    nx, ny = x + dx, y + dy

                    # Дополнительная проверка на занятость клетки
                    if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                        new_grid[nx, ny] = 2  # Перемещаем антибиотик в новую клетку
                        new_grid[x, y] = 0  # Обнуляем старую клетку

            elif grid[x, y] == 3:  # Если в клетке есть антибиотик 3
                possible_moves = []  # Список возможных направлений для движения антибиотика

                # Проверка соседей (вверх, вниз, влево, вправо)
                if x - 1 >= 0 and grid[x - 1, y] == 0:
                    possible_moves.append((-1, 0))  # Вверх
                if x + 1 < grid_size and grid[x + 1, y] == 0:
                    possible_moves.append((1, 0))  # Вниз
                if y - 1 >= 0 and grid[x, y - 1] == 0:
                    possible_moves.append((0, -1))  # Влево
                if y + 1 < grid_size and grid[x, y + 1] == 0:
                    possible_moves.append((0, 1))  # Вправо

                # Если есть возможные движения, случайным образом выбрать одно
                if possible_moves:
                    dx, dy = random.choice(possible_moves)
                    nx, ny = x + dx, y + dy

                    # Дополнительная проверка на занятость клетки
                    if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                        new_grid[nx, ny] = 3  # Перемещаем антибиотик в новую клетку
                        new_grid[x, y] = 0  # Обнуляем старую клетку

    # Обновляем сетку после выполнения всех шагов диффузии
    grid = new_grid

    # Обновляем визуализацию
    im.set_data(grid)
    plt.draw()  # Отрисовываем обновление
    plt.pause(0.1)  # Пауза для обновления

# Визуализация процесса с блокировкой
for frame in range(100):
    update(frame)

plt.show()
