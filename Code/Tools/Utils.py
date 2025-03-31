def find_closest(data:list, value:float):
    return min(data, key=lambda x:abs(x-float(value)))