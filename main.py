#Created by João Gabriel Barbosa de Luna - iLIKA/UFPE

#import omicscope as omics
import os
from datetime import datetime


def create_dir(directory, project_name, username):
    os.chdir(directory)
    os.mkdir(project_name)
    os.mkdir(os.path.join(project_name, "tables"))
    os.mkdir(os.path.join(project_name, "plots"))
    with open(os.path.join(project_name, f"{project_name}.log"), "w+") as log_file:
        log_file.write(f"Project name: {project_name}\n")
        log_file.write(f"Owned by: {username}\n")
        log_file.write(f"Created on: {datetime.now()}")

def read_proteomics_file(filepath, method, control):
    if not os.path.exists(filepath):
        raise FileNotFoundError("Arquivo não encontrado. Verifique o caminho inserido e tente novamente.")
    elif os.path.getsize(filepath) == 0:
        raise Exception("O arquivo fornecido está vazio. Verifique o conteúdo do mesmo e tente novamente.")
    else:
        data = omics.OmicScope(filepath, Method=method, ControlGroup=control)
    return data


create_dir("C:/Users/gabri/Documents/proteomics", "test", "joão")



