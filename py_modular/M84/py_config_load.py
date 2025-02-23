import re
import ast


def load_configuration(conf_file):
    conf_var = {}
    with open(conf_file, "r") as file:
        for line in file:
            line = line.split("#")[0].strip()  # Eliminar comentarios y espacios
            if line:  # Ignorar líneas vacías
                match = re.match(r"^([A-Z_]+)\s*=\s*(.+)", line)  # Buscar VARIABLE = VALOR
                if match:
                    key, value = match.groups()
                    value = value.strip()
                    
                    # Convertir valores: int si es número, lista o tupla si tiene formato correcto
                    try:
                        value = ast.literal_eval(value)  # Convierte '0' a 0, '["x"]' a lista, etc.
                    except (ValueError, SyntaxError):
                        pass  # Si no es un número o estructura reconocible, se deja como string
                    
                    conf_var[key] = value
    
    return conf_var