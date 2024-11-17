import matplotlib.pyplot as plt
import numpy as np
from random import randint
from Bio import SeqIO
import Levenshtein

# Укажите путь к вашему файлу FNA/FASTA
file_path = r"C:\Users\Андрей\Desktop\GCA_000008865.2_ASM886v2_genomic.fna"

# Чтение референсного генома
with open(file_path, "r") as file:
    record = next(SeqIO.parse(file, "fasta"))  # Считываем первую запись из файла
    ref_genome = str(record.seq)  # Преобразуем последовательность в строку

nucleotides = ['A', 'T', 'C', 'G']
class Cell:

    dna = ''

    is_alive = True
    HP = 100

    def __init__(self, dna:str, coords_x:int, coords_y:int, conc_antibiotics:int):
        self.dna = dna
        self.X = coords_x
        self.Y = coords_y
        self.is_alive = True
        self.conc_antibiotics = conc_antibiotics
        self.age_epoch = 0

    def mutate(self, num_muts):
        dna_list = list(self.dna)  # Преобразуем строку ДНК в список для удобства изменений
        for i in range(num_muts):
            position = randint(0, len(dna_list) - 1)  # Выбираем случайную позицию в ДНК
            # Выбираем случайный нуклеотид, индексируем nucleotides от 0 до 3
            dna_list[position] = nucleotides[randint(0, len(nucleotides) - 1)]
        return dna_list # Возвращаем мутированную ДНК в виде строки

    def divide(self):
        #print('2: ', end = '')
        #start = time.time()
        first_dna = self.mutate(1)
        second_dna = self.mutate(1)
        return first_dna, second_dna
        #print(time.time() - start)
        # сделать чтобы клетки делились, график изменения количества клеток
        # бахнуть график, показывающий изменение гетерогенности генофонда популяции

    def check_alive(self, resist_gen):
        """
        Проверка, жива ли клетка, в зависимости от концентрации антибиотика и наличия гена резистентности.
        """
        k = 0.1
        prob_resist = np.exp(-k*self.conc_antibiotics)
        if resist_gen in self.dna:
            #print("Cell is alive")
            self.age_epoch += 1
            self.alive = True
        elif resist_gen not in self.dna and prob_resist > 0.7:
            #print("Cell is alive")
            self.age_epoch += 1
            self.alive = True
        elif prob_resist <=0.7:
            #print("Cell is dead")
            self.alive = False
        return self.alive, self.age_epoch


class Grid_cell:

    is_celled = False
    X = -1
    Y = -1

    def __init__(self, coords_x, coords_y):
        self.X = coords_x
        self.Y = coords_y

    # Группа 1: Методы, связанные с состоянием клетки (жива ли она и её смерть)

    # def fragment_dna(self):
    #     """
    #     Делит ДНК клетки на случайные фрагменты, если клетка погибла.
    #     """
    #     dna_length = len(self.dna)
    #     fragments = []
    #     i = 0
    #     while i < dna_length:
    #         frag_length = np.random.randint(1, 5)  # Длина фрагмента от 1 до 5 нуклеотидов
    #         fragments.append(self.dna[i:i + frag_length])
    #         i += frag_length
    #     print("ДНК фрагментирована:", fragments)
    #     return fragments

    # # Группа 2: Методы, связанные с дубликацией и трансформацией ДНК
    # def duplication(self, num_of_muts):
    #     temp_res = list(self.dna)
    #     position = np.random.randint(0, len(temp_res)-1)
    #     temp_res[position] = choice(['A', 'T', 'G', 'C'], num_of_muts)
    #     dupl_dna = temp_res
    #     return temp_res, dupl_dna


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


    # def diffuse_antibiotics(self, grid, diffusion_rate=0.1):
    #     """
    #     Моделирование диффузии антибиотиков на сетке.
    #     :param grid: Двумерный массив (numpy.ndarray), представляющий концентрации антибиотиков.
    #     :param diffusion_rate: Коэффициент диффузии (доля перемещаемого вещества).
    #     :return: Обновленная сетка.
    #     """
    #     new_grid = grid.copy()
    #     grid_size = grid.shape[0]

    #     for x in range(grid_size):
    #         for y in range(grid_size):
    #             if grid[x, y] > 0:  # Если в ячейке есть антибиотик
    #                 diffusion_amount = grid[x, y] * diffusion_rate

    #                 # Распределяем часть концентрации на соседей (4 соседние ячейки)
    #                 for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
    #                     nx, ny = x + dx, y + dy
    #                     if 0 <= nx < grid_size and 0 <= ny < grid_size:
    #                         new_grid[nx, ny] += diffusion_amount / 4

    #                 # Уменьшаем текущую концентрацию
    #                 new_grid[x, y] -= diffusion_amount

    #     return new_grid


# Размер сетки
grid_size = 100  # Размер сетки 10x10
grid = np.zeros((grid_size, grid_size))  # Инициализация пустой сетки (0 - пустое место)

# Количество клеток и антибиотиков
num_organisms = 3000  # Количество клеток
num_antibiotics = 7000  # Количество антибиотиков

# Создание клеток
cells = []  # Список для хранения клеток
antibiotics = []  # Список для хранения клеток
heterogenity = []
epoch = []

def new_antibiotic(x,y,grid):
    antibiotic = Grid_cell(x, y)  # Создаем клетку
    antibiotics.append(antibiotic)  # Добавляем клетку в список
    grid[x][y] = 2  # Размещение клетки на сетке
    return grid[x][y]

# Генерация случайных позиций для клеток (1)
for _ in range(num_antibiotics):
    x, y = np.random.randint(0, grid_size, 2)  # Случайные координаты x и y для клетки
    new_antibiotic(x,y,grid)

def new_organism(x,y,ref_genome,grid,num_muts):
    count_antibiotic = 0
    # Проверяем 8 соседних клеток вокруг клетки
    for dx in range(-1, 2):  # x-сдвиг от -1 до 1
        for dy in range(-1, 2):  # y-сдвиг от -1 до 1
            # Пропускаем саму клетку
            if dx == 0 and dy == 0:
                continue
            # Проверяем, чтобы не выйти за границы
            if 0 <= x + dx < grid_size and 0 <= y + dy < grid_size:
                # Если соседняя клетка содержит антибиотик (2)
                if grid[x + dx][y + dy] == 2:
                    count_antibiotic += 1

    org_dna = ''.join(np.random.choice(nucleotides, size=10))  # Пример ДНК длиной 10
    organism = Cell(org_dna, x, y, count_antibiotic)
    organism.dna = organism.mutate(num_muts)
    cells.append(organism)  # Добавляем клетку в список
    grid[x][y] = 1  # Размещение клетки на сетке
    return grid[x][y]

# Генерация случайных позиций для клеток (1)
for _ in range(num_organisms):
    x, y = np.random.randint(0, grid_size, 2)  # Случайные координаты x и y для клетки
    dna = ''.join(np.random.choice(nucleotides, size=10))
    new_organism(x,y,dna,grid,5)

start = np.shape(cells)


def chek_neighbours_grid(x,y):
    possible_moves = []
    # Проверка соседей (вверх, вниз, влево, вправо)
    if x - 1 >= 0 and grid[x - 1, y] == 0:
        possible_moves.append((-1, 0))  # Вверх
    if x + 1 < grid_size and grid[x + 1, y] == 0:
        possible_moves.append((1, 0))  # Вниз
    if y - 1 >= 0 and grid[x, y - 1] == 0:
        possible_moves.append((0, -1))  # Влево
    if y + 1 < grid_size and grid[x, y + 1] == 0:
        possible_moves.append((0, 1))  # Вправо
    return possible_moves

def chek_resist(resist_gen, new_grid):
    gen = resist_gen  # Пример гена резистентности
    for cell in cells:
        is_alive, is_age_epoch = cell.check_alive(gen)
        if not is_alive:
            new_grid[cell.X][cell.Y] = 0
            cells.remove(cell)
        if is_age_epoch > 10:
            print("10 EPOCH died")
            new_grid[cell.X][cell.Y] = 0
            cells.remove(cell)

def calculate_heterogeneity(cells):
    """
    Рассчитывает гетерогенность генофонда.

    Args:
        cells (list): Список клеток.

    Returns:
        float: Гетерогенность генофонда G.
    """
    # Получаем только DNA последовательности клеток
    genomes = [cell.dna for cell in cells]

    n = len(genomes) - 1  # Количество пар геномов
    if n <= 0:
        return 0  # Если меньше двух геномов, гетерогенность = 0

    total_distance = 0

    for i in range(n):
        g_i = Levenshtein.distance(genomes[i], genomes[i+1])
        total_distance += g_i / len(genomes[i])  # Нормализация расстояния

    G = total_distance / n  # Усреднение
    return G


# Настройка визуализации
fig, ax = plt.subplots(figsize=(8, 8))
cmap = plt.get_cmap("viridis")
ax.set_title("Карта клеток и антибиотиков")
im = ax.imshow(grid, cmap=cmap, interpolation='nearest')

def update(EPOCH):
    global grid
    new_grid = grid.copy()  # Копируем текущую сетку для обновлений

    if EPOCH %5 !=0:
        # Один шаг диффузии для антибиотиков и клеток
        for x in range(grid_size):
            for y in range(grid_size):
                if grid[x, y] == 1:  # Если в клетке есть клетка (cell)
                    possible_moves = chek_neighbours_grid(x,y)
                    # Если есть возможные движения, случайным образом выбрать одно
                    if possible_moves:
                        dx, dy = possible_moves[np.random.choice(len(possible_moves))]  # Правильный выбор
                        nx, ny = x + dx, y + dy
                        # Дополнительная проверка на занятость клетки
                        if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                            for cell in cells:
                                if cell.X == x and cell.Y == y:
                                    new_organism(nx,ny,cell.dna,new_grid,0)
                                    cells.remove(cell)
                                    new_grid[x][y] = 0

                if grid[x, y] == 2:  # Если в клетке есть клетка (antibiotic)
                    possible_moves = chek_neighbours_grid(x,y)
                    # Если есть возможные движения, случайным образом выбрать одно
                    if possible_moves:
                        dx, dy = possible_moves[np.random.choice(len(possible_moves))]  # Правильный выбор
                        nx, ny = x + dx, y + dy

                        # Дополнительная проверка на занятость клетки
                        if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                            new_antibiotic(nx,ny,new_grid)  # Добавляем клетку в список
                            new_grid[x][y] = 0  # Перемещаем клетку в новую клетку
    chek_resist("ATGC", new_grid)

    if EPOCH%5 == 0:
         for x in range(grid_size):
             for y in range(grid_size):
                 if grid[x, y] == 1:  # Если в клетке есть клетка (cell)
                     possible_moves = chek_neighbours_grid(x,y)
                     # Если есть возможные движения, случайным образом выбрать одно
                     if possible_moves:
                         dx, dy = possible_moves[np.random.choice(len(possible_moves))]  # Правильный выбор
                         nx, ny = x + dx, y + dy
                         # Дополнительная проверка на занятость клетки
                         if new_grid[nx, ny] == 0:  # Проверка, что клетка пуста
                             for cell in cells:
                                 if cell.X == x and cell.Y == y:
                                     new_organism(nx,ny,cell.mutate(1),new_grid,0)
    chek_resist("ATGC", new_grid)


    # Обновляем сетку после выполнения всех шагов диффузии
    grid = new_grid
    im.set_data(grid)
    plt.draw()  # Отрисовываем обновление
    plt.pause(0.1)  # Пауза для обновления

    for cell in cells:
        n = len(cells)/2
        x_i = len(cell.dna)

    G = calculate_heterogeneity(cells)
    heterogenity.append(G)
    epoch.append(EPOCH)

# Визуализация процесса с блокировкой
for EPOCH in range(50):
    update(EPOCH)

plt.show()
print(start)
print(np.shape(cells))


plt.plot(epoch,heterogenity)
plt.legend("Heterogenity")
plt.show()
