#Created by João Gabriel Barbosa de Luna - iLIKA/UFPE
import json
import omicscope as omics
import pandas as pd
import os
from datetime import datetime
import subprocess


#---------------------------------------------------------------------------------#
def create_dir(targetdirectory, projectname, username):
    os.chdir(targetdirectory)
    os.mkdir(projectname)
    os.mkdir(os.path.join(projectname, "tables"))
    os.mkdir(os.path.join(projectname, "plots"))
    with open(os.path.join(projectname, f"{projectname}.log"), "w+") as logfile:
        logfile.write(f"Project name: {projectname}\n")
        logfile.write(f"Owned by: {username}\n")
        logfile.write(f"Created on: {datetime.now()}")
    return os.path.join(targetdirectory, projectname)

def read_proteomics_file(rawfilepath, method, control):
    if not os.path.exists(rawfilepath):
        raise FileNotFoundError("Arquivo não encontrado. Verifique o caminho inserido e tente novamente.")
    elif os.path.getsize(rawfilepath) == 0:
        raise Exception("O arquivo fornecido está vazio. Verifique o conteúdo do mesmo e tente novamente.")
    else:
        data = omics.OmicScope(rawfilepath, Method=method, ControlGroup=control)
    return data

def create_xlsx_file(dataframe, name, tablesdir):
    filename = os.path.join(tablesdir, f"{name}.xlsx")
    dataframe.to_excel(filename)

def create_json_file(dictionary, name, tablesdir):
    filename = os.path.join(tablesdir, f"{name}.json")
    with open(filename, "w+") as js:
        json.dump(dictionary, js)

#----------------------------------------------------------------------------------#
root = os.path.dirname(os.path.abspath(__file__))
#----------------------------------------------------------------------------------#

if __name__ == "__main__":
    while True:
        try:
            print("Proteoanalyzer 1.0")
            # Gerando diretórios e arquivo de log das análises
            user_name = input("Insira o nome do usuário: ")
            project_name = input("Insira um nome para o projeto a ser analisado: ")
            target_directory = input("Insira o caminho onde os arquivos da análise serão criados (ex: C:/Users/user/Documents): ")
            directory = create_dir(target_directory, project_name, user_name)
            log_file = os.path.join(directory, f"{project_name}.log")
            tables_dir = os.path.join(directory, "tables")
            plots_dir = os.path.join(directory, "plots")
            subprocess.run(["python", "tables.py", tables_dir, plots_dir], cwd=root)

            #Obtendo dados das análises a partir do usuario
            raw_file_path = input("Insira o caminho da planilha com os dados da proteômica (não esqueça de incluir a extensão do arquivo)")
            proteomics_method = input("Insira o software de análise dos dados brutos (MaxQuant, Progenesis, General, DIA-NN, PatternLab)")
            control_group = input("Insira o nome do grupo controle, exatamente como está na planilha: ")

            #Lendo os dados e criando os arquivos com as planilhas
            raw_file_data = read_proteomics_file(raw_file_path, proteomics_method, control_group)
            params_dataframe = pd.Dataframe(raw_file_data.Params())
            create_xlsx_file(params_dataframe, "Parâmetros", tables_dir)
            conditions_dictionary = {"Condições": f"{raw_file_data.Conditions}",
                                     "Controle": f"{raw_file_data.ControlGroup}"}
            create_json_file(conditions_dictionary, "Condições do estudo", tables_dir)
            all_data_dataframe = pd.DataFrame(raw_file_data.quant_data())
            create_xlsx_file(all_data_dataframe, "Dados brutos", tables_dir)
            deps_dataframe = pd.DataFrame(raw_file_data.deps())
            create_xlsx_file(deps_dataframe, "DEPs", tables_dir)

            #Plotando os gráficos
            raw_file_data.bar_ident(logscale=True,
                                    dpi=300,
                                    save=os.path.join(plots_dir, f"{project_name}_Barplot_Identificação.tiff"))

            dynamic_range_selection = input("Deseja plotar o gráfico de Dynamic range? (y/n): ").lower()
            while True:
                if dynamic_range_selection == "y":
                    dynamic_range_proteins = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula e os nomes entre aspas simples: ")
                    dynamic_range_proteins_list = [p.strip() for p in dynamic_range_proteins.strip(",")]
                    raw_file_data.DynamicRange(*dynamic_range_proteins_list,
                                               dpi=300,
                                               save=os.path.join(plots_dir, f"{project_name}_Dynamic_range.tiff"))
                    break
                elif dynamic_range_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {dynamic_range_selection}")
                    continue

            #Volcano_plot_selection = input("Deseja plotar o volcano plot? (y/n): ")

            MA_plot_selection = input("Deseja plotar o grático MA plot? (y/n): ")
            while True:
                if MA_plot_selection == "y":
                    MA_plot_proteins = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula e os nomes entre aspas simples: ")
                    MA_plot_proteins_list = [p.strip() for p in MA_plot_proteins.strip(",")]
                    raw_file_data.MAplot(*MA_plot_proteins_list,
                                         dpi=300,
                                         save=os.path.join(plots_dir, f"{project_name}_MA_plot.tiff"))
                    break
                elif MA_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {MA_plot_selection}")
                    continue

            Normalization_plot_selection = input("Deseja plotar o gráfico de normalização? (y/n): ")
            while True:
                if Normalization_plot_selection == "y":
                    raw_file_data.normalization_boxplot(dpi=300,
                                                        save=os.path.join(plots_dir, f"{project_name}_Normalização_boxplot.tiff"))
                    break
                elif Normalization_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {Normalization_plot_selection}")
                    continue

            Proteins_bar_plot_selection = input("Deseja plotar o gráfico de abundância de proteinas entre as condições do estudo? (y/n): ")
            while True:
                if Proteins_bar_plot_selection == "y":
                    Proteins_bar_plot_proteins = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula e os nomes entre aspas simples: ")
                    Proteins_bar_plot_proteins_list = [p.strip() for p in Proteins_bar_plot_proteins.strip(",")]
                    raw_file_data.bar_protein(*Proteins_bar_plot_proteins_list,
                                         palette='viridis',
                                         dpi=300,
                                         save=os.path.join(plots_dir, f"{project_name}_Bar_Plot_Proteínas_Entre_Condições.tiff"))
                    break
                elif Proteins_bar_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {Proteins_bar_plot_selection}")
                    continue

            Proteins_boxplot_selection = input("Deseja plotar o gráfico boxplot de abundância de proteinas entre as condições do estudo? (y/n): ")
            while True:
                if Proteins_boxplot_selection == "y":
                    Proteins_boxplot_proteins = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula e os nomes entre aspas simples: ")
                    Proteins_boxplot_proteins_list = [p.strip() for p in Proteins_boxplot_proteins.strip(",")]
                    raw_file_data.boxplot_protein(*Proteins_boxplot_proteins_list,
                                              palette='viridis',
                                              dpi=300,
                                              save=os.path.join(plots_dir, f"{project_name}_Boxplot_Proteínas_Entre_Condições.tiff"))
                    break
                elif Proteins_boxplot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {Proteins_boxplot_selection}")
                    continue

            Heatmap_selection = input("Deseja plotar o heatmap do nível de expressão das proteínas? (y/n): ")
            while True:
                if Heatmap_selection == "y":
                    raw_file_data.heatmap(dpi=300,
                                          linewidth=0,
                                          save=os.path.join(plots_dir, f"{project_name}_Heatmap_expressão.tiff"))
                    break
                elif Heatmap_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {Heatmap_selection}")
                    continue

            Correlation_plot_selection = input("Deseja plotar o heatmap de correlação entre os grupos? (y/n): ")
            while True:
                if Correlation_plot_selection == "y":
                    raw_file_data.correlation(dpi=300,
                                          linewidth=0,
                                          save=os.path.join(plots_dir, f"{project_name}_Heatmap_correlação.tiff"))
                    break
                elif Correlation_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {Correlation_plot_selection}")
                    continue

            PCA_plot_selection = input("Deseja plotar o gráfico de PCA? (y/n): ")
            while True:
                if PCA_plot_selection == "y":
                    raw_file_data.pca(pvalue=0.05,
                                      dpi=300,
                                      save=os.path.join(plots_dir, f"{project_name}_PCA.tiff"))
                    break
                elif PCA_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {PCA_plot_selection}")
                    continue

            K_mean_plot_selection = input("Deseja plotar o gráfico de K-means? (y/n): ")
            while True:
                if K_mean_plot_selection == "y":
                    raw_file_data.k_trend(dpi=300,
                                          save=os.path.join(plots_dir, f"{project_name}_Kmeans.tiff"))
                    break
                elif K_mean_plot_selection == "n":
                    break
                else:
                    print(f"Valor inválido inserido: {K_mean_plot_selection}")
                    continue

        except ValueError as error:
            print(f"O valor inserido é invalido: {error}")
        except Exception as error:
            print(f"Erro inesperado: {error}")








