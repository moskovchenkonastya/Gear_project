# Gear_project

Задача: разработать программу, которая находит и заменяет сломанные шестеренки на изображении механизма

1) Бинаризация изображения, с помощью отсечения пикселей по порогу
2) Выделение связных компонент на изображении
3) Классификация связных компонент
4) Вычисление геометрических признаков связных компонент
- положение центра окружности 
- радиус меньшей окружности, образующей шестеренку
- радиус большей окружности, образующей шестеренку
5) Выбор шестеренки, которую необходимо установить на ось
6) Поиск сломанной шестеренки с помощью алгоритма distance transform
7) Определения колличества зубцов на каждой шестеренке

# Cтруктура программы
