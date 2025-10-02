import numpy as np
from scipy.optimize import root
from math import degrees

def z(x, Rc, K):
    """
    Функция, описывающая кривую второго порядка (коническую поверхность).
    """
    return (x**2) / (Rc * (1 + np.sqrt(1 - (1 + K) * (x**2 / Rc**2))))

def equations(vars, D, ZR, Rc, K):
    """
    Система уравнений:
    1. (x1 - x2)**2 + (z(x1) - z(x2))**2 = D**2
    2. Расстояние между параллельными прямыми = ZR
       Одна прямая: проходит через центр отрезка и перпендикулярна ему
       Вторая прямая: проходит через (0,0) и параллельна первой
    """
    x1, x2 = vars
    
    # Вычисляем значения функции z(x) в точках x1 и x2
    z1 = z(x1, Rc, K)
    z2 = z(x2, Rc, K)
    
    # Координаты центра отрезка
    x_center = (x1 + x2) / 2
    z_center = (z1 + z2) / 2
    
    # Вектор направления отрезка (от x1,z1 к x2,z2)
    dx = x2 - x1
    dz = z2 - z1
    
    nx = dx  # Перпендикуляр к направлению отрезка
    nz = dz
    
    # Уравнение прямой, проходящей через центр и перпендикулярной отрезку:
    # (x - x_center) * nx + (z - z_center) * nz = 0
    # => nx*x + nz*z - (nx*x_center + nz*z_center) = 0
    
    # Расстояние от точки (0,0) до этой прямой:
    # distance = |A*0 + B*0 + C| / sqrt(A^2 + B^2)
    # где A = nx, B = nz, C = -(nx*x_center + nz*z_center)
    A = nx
    B = nz  
    C = -(nx * x_center + nz * z_center)
    
    # Это расстояние должно равняться ZR
    calculated_distance = abs(C) / np.sqrt(A**2 + B**2)
    
    # Первое уравнение: расстояние между точками = D
    eq1 = (x1 - x2)**2 + (z1 - z2)**2 - D**2
    
    # Второе уравнение: расстояние между параллельными прямыми = ZR
    eq2 = calculated_distance - ZR
    
    return [eq1, eq2]

def solve_system(D, ZR, Rc, K, method='lm'):
    """
    Решение системы уравнений.
    """
    # Начальное приближение
    initial_guess = [ZR - D/2, ZR + D/2]
    
    # Решаем систему уравнений
    solution = root(equations, initial_guess, args=(D, ZR, Rc, K), method=method)
    
    if solution.success:
        return solution.x
    else:
        print(f"Решение не найдено. Сообщение: {solution.message}")
        return None

def equations_center(vars, Rc, K, A, B, C):
    """
    Система уравнений:
    1. z = f(x) уравнение кривой второго порядка
    2. z = g(x) уравнение оси внеосевой апертуры
    """
    xx, zz = vars
    
    # Вычисляем невязку zz - z(x) в точке xx
    eq1 = zz - z(xx, Rc, K)
    eq2 = A*xx + B*zz + C  # уравнение прямой    
    
    return [eq1, eq2]

def solve_center(D, ZR, Rc, K, method='lm'):
    """
    Решение системы уравнений.
    """
    x1, x2 = solve_system(D, ZR, Rc, K, method)
    z1, z2 = z(x1, Rc, K), z(x2, Rc, K)

    # Координаты центра отрезка
    x_center = (x1 + x2) / 2
    z_center = (z1 + z2) / 2

    # уравнение прямой, проходящей через центр и перпендикулярной отрезку:
    A = x2 - x1
    B = z2 - z1
    C = -(A * x_center + B * z_center)

    # Начальное приближение
    initial_guess = [ZR, z(ZR, Rc, K)]
    
    # Решаем систему уравнений
    solution = root(equations_center, initial_guess, args=(Rc, K, A, B, C), method=method)
    
    if solution.success:
        return solution.x
    else:
        print(f"Решение не найдено. Сообщение: {solution.message}")
        return None

def radians_to_dms(radians, decimals=2):
    """
    Преобразует угол из радианов в строку формата "градусы° минуты' секунды""
    
    Параметры:
        radians : float
            Угол в радианах
        decimals : int, optional
            Количество знаков после запятой для секунд (по умолчанию 2)
    
    Возвращает:
        str
            Строка в формате "X° Y' Z.ZZ""
    """
    # Переводим радианы в градусы
    degrees_total = degrees(radians)
    
    # Определяем знак
    sign = -1 if degrees_total < 0 else 1
    degrees_total = abs(degrees_total)
    
    # Выделяем целые градусы
    degrees_int = int(degrees_total)
    
    # Вычисляем минуты
    minutes_total = (degrees_total - degrees_int) * 60
    minutes_int = int(minutes_total)
    
    # Вычисляем секунды
    seconds = (minutes_total - minutes_int) * 60
    
    # Округляем секунды до указанного количества знаков
    seconds = round(seconds, decimals)
    
    # Обрабатываем перенос разрядов (если секунды = 60)
    if seconds >= 60:
        seconds = 0
        minutes_int += 1
        if minutes_int >= 60:
            minutes_int = 0
            degrees_int += 1
    
    # Форматируем секунды
    if decimals == 0:
        seconds_str = f"{int(seconds)}"
    else:
        seconds_str = f"{seconds:0.{decimals}f}"
    
    # Формируем итоговую строку со знаком
    sign_str = "-" if sign == -1 else ""
    
    return f"{sign_str}{degrees_int}° {minutes_int}' {seconds_str}\""

# Пример использования с вашими данными
if __name__ == "__main__":
    Rc = 1294.2
    K = -0.23455
    D = 214.0
    ZR = 166.0
    
    result = solve_system(D, ZR, Rc, K)
    center = solve_center(D, ZR, Rc, K)

    if result is not None:
        x1_sol, x2_sol = result
        XC, ZC = center if center is not None else (None, None)

        print(f"Решение найдено:")
        print(f"x1 = {x1_sol:.6f}")
        print(f"x2 = {x2_sol:.6f}")
        
        # Проверка решения
        z1 = z(x1_sol, Rc, K)
        z2 = z(x2_sol, Rc, K)
        
        # Проверка первого уравнения
        actual_distance = np.sqrt((x1_sol - x2_sol)**2 + (z1 - z2)**2)
        
        # Проверка второго уравнения
        x_center = (x1_sol + x2_sol) / 2
        z_center = (z1 + z2) / 2
        dx = x2_sol - x1_sol
        dz = z2 - z1
        length = np.sqrt(dx**2 + dz**2)
        nx = dx / length if length > 0 else 1
        nz = dz / length if length > 0 else 0
        A = nx
        B = nz
        C = -(nx * x_center + nz * z_center)
        calculated_distance = abs(C) / np.sqrt(A**2 + B**2)
        dxc, dzc = XC - x_center, ZC - z_center
        length_c = np.sqrt(dxc**2 + dzc**2)
        dx_c, dz_c = dxc / length_c, dzc / length_c

        print(f"\nПроверка решений:")
        print(f"Заданное расстояние D = {D:.6f}")
        print(f"Фактическое расстояние между точками = {actual_distance:.6f}")
        print(f"Заданное расстояние между прямыми ZR = {ZR:.6f}")
        print(f"Фактическое расстояние между прямыми = {calculated_distance:.6f}")
        
        print(f"\nДополнительная информация:")
        print(f"z(x1) = {z1:.6f}")
        print(f"z(x2) = {z2:.6f}")
        print(f"Координаты центра отрезка: ({x_center:.6f}, {z_center:.6f})")
        print(f"Координаты центра апертуры: XC = {XC:.6f}, ZC = {ZC:.6f}\nПроверка решения: z(XC) = {z(XC, Rc, K):.6f}")
        print(f"Вектор нормали: ({dx_c:.6f}, {dz_c:.6f})")
        print(f"Вектор направления: ({dx/length:.6f}, {dz/length:.6f})")
        print(f"Ортогональность векторов: {dx_c*dx + dz_c*dz:.6f} (должно быть ~0)")
        print(f"Угол наклона: {radians_to_dms(np.arctan(dz/dx))}")
    else:
        print("Решение не найдено. Попробуйте другие параметры или метод решения.")