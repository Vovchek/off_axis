import numpy as np
from scipy.optimize import root
import math

# Определение функции z(x) для кривой второго порядка
def z(x, Rc, K):
    """
    Функция, описывающая кривую второго порядка (коническую поверхность).
    
    Параметры:
        x : array_like
            Координата x (может быть скаляром или массивом)
        Rc : float
            Радиус вершинной кривизны
        K : float
            Коническая постоянная
    
    Возвращает:
        z : array_like
            Значение функции z в точке x
    """
    return (x**2) / (Rc * (1 + np.sqrt(1 - (1 + K) * (x**2 / Rc**2))))

# Система уравнений для решения
def equations(vars, D, ZR, Rc, K):
    """
    Система уравнений:
    (x1 - x2)**2 + (z(x1) - z(x2))**2 = D**2
    (x1 + x2)**2 + (z(x1) + z(x2))**2 = 4 * ZR
    
    Параметры:
        vars : list
            Вектор переменных [x1, x2]
        D : float
            Параметр D из уравнений
        ZR : float
            Параметр ZR из уравнений
        Rc : float
            Радиус вершинной кривизны для функции z(x)
        K : float
            Коническая постоянная для функции z(x)
    
    Возвращает:
        list
            Вектор невязок уравнений
    """
    x1, x2 = vars
    
    # Вычисляем значения функции z(x) в точках x1 и x2
    z1 = z(x1, Rc, K)
    z2 = z(x2, Rc, K)
    
    # Первое уравнение: (x1 - x2)^2 + (z1 - z2)^2 - D^2 = 0
    eq1 = (x1 - x2)**2 + (z1 - z2)**2 - D**2
    
    # Второе уравнение: (x1 + x2)^2 + (z1 + z2)^2 - 4*ZR = 0
    eq2 = (x1 + x2)**2 + (z1 + z2)**2 - 4 * ZR * ZR
    
    return [eq1, eq2]

# Функция для решения системы с интеллектуальным начальным приближением
def solve_system(D, ZR, Rc, K, method='lm'):
    """
    Решение системы уравнений с логичным начальным приближением.
    
    Параметры:
        D : float
            Параметр D (расстояние между точками)
        ZR : float
            Параметр ZR (связан с положением центра)
        Rc : float
            Радиус вершинной кривизны
        K : float
            Коническая постоянная
        method : str
            Метод решения ('lm', 'hybr', etc)
    
    Возвращает:
        tuple (x1, x2) или None в случае неудачи
    """
    # Логичное начальное приближение: точки симметричны относительно ZR
    initial_guess = [ZR - D/2, ZR + D/2]
    
    # Решаем систему уравнений
    solution = root(equations, initial_guess, args=(D, ZR, Rc, K), method=method)
    
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
    degrees_total = math.degrees(radians)
    
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

# Пример использования
if __name__ == "__main__":
    # Задаем параметры
    Rc = 1294.2  # Радиус вершинной кривизны
    K = -0.23455    # Коническая постоянная
    D = 214.0    # Параметр D (расстояние)
    ZR = 166.0   # Параметр ZR (положение центра)
    
    # Решаем систему
    result = solve_system(D, ZR, Rc, K)
    
    if result is not None:
        x1_sol, x2_sol = result
        print(f"Решение найдено:")
        print(f"x1 = {x1_sol:.6f}")
        print(f"x2 = {x2_sol:.6f}")
        
        # Проверка решения
        z1 = z(x1_sol, Rc, K)
        z2 = z(x2_sol, Rc, K)
        
        # Вычисляем фактические расстояния
        actual_distance = np.sqrt((x1_sol - x2_sol)**2 + (z1 - z2)**2)
        actual_center_distance = np.sqrt(((x1_sol + x2_sol)/2)**2 + ((z1 + z2)/2)**2)
        
        print(f"\nПроверка решений:")
        print(f"Заданное расстояние D = {D:.6f}")
        print(f"Фактическое расстояние между точками = {actual_distance:.6f}")
        print(f"Заданное расстояние до центра ZR = {ZR:.6f}")
        print(f"Фактическое расстояние до центра = {actual_center_distance:.6f}")
        
        print(f"\nНевязки уравнений:")
        eq1_residual = (x1_sol - x2_sol)**2 + (z1 - z2)**2 - D**2
        eq2_residual = (x1_sol + x2_sol)**2 + (z1 + z2)**2 - 4 * ZR**2
        print(f"Уравнение 1: {eq1_residual:.2e}")
        print(f"Уравнение 2: {eq2_residual:.2e}")
        
        # Дополнительная информация
        print(f"\nДополнительная информация:")
        print(f"z(x1) = {z1:.6f}")
        print(f"z(x2) = {z2:.6f}")
        print(f"Координаты центра: ({ (x1_sol + x2_sol)/2:.6f}, { (z1 + z2)/2:.6f})")
        print(f"Угол наклона: {radians_to_dms(np.arctan((z2-z1)/(x2_sol-x1_sol)))}")
    else:
        print("Решение не найдено. Попробуйте другие параметры или метод решения.")