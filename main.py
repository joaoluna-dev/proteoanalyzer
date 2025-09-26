# Created by João Gabriel Barbosa de Luna - iLIKA/UFPE
import json
import omicscope as omics
import pandas as pd
import os
from datetime import datetime
import math
import logging
import shutil


# ---------------------------------------------------------------------------------#
def create_dir(targetdirectory, projectname):
    """Cria a estrutura de diretórios do projeto."""
    if not os.path.exists(targetdirectory):
        logging.error(f"O diretório {targetdirectory} inserido não existe.")
        raise FileNotFoundError(f"O diretório de destino não existe: {targetdirectory}")

    project_path = os.path.join(targetdirectory, projectname)

    # Verifica se o diretório do projeto já existe
    if os.path.exists(project_path):
        response = input(f"O diretório '{projectname}' já existe no caminho selecionado. Deseja sobrescrever? (y/n): ").lower().strip()
        if response == 'n':
            while True:
                another_dir = input("Deseja selecionar outro diretório? (y/n): ").lower().strip()
                if another_dir == "y":
                    another_directory = input("Insira o caminho onde os arquivos da análise serão criados (ex: C:/Users/user/Documents): ").strip().strip('"')
                    project_path = os.path.join(another_directory, projectname)
                    logging.info(f"Novo diretório criado: {project_path}")
                    if os.path.exists(project_path):
                        print(f"O diretório {projectname} também existe no caminho selecionado.")
                        logging.error(f"O diretório {projectname} existe no caminho selecionado.")
                        continue
                    else:
                        logging.info(f"O diretório {projectname} foi criado em {another_directory}")
                        print(f"O diretório {projectname} foi criado em {another_directory}")
                        break
                elif another_dir == "n":
                    logging.info("Operação cancelada pelo usuário.")
                    raise Exception("Operação cancelada pelo usuário.")
                else:
                    print(f"Opção inválida inserida: {another_dir}")
                    continue
        elif response == 'y':
            shutil.rmtree(project_path)
            os.mkdir(project_path)
            logging.info(f"O diretório {projectname} foi recriado em {targetdirectory}")
            print(f"O diretório {projectname} foi recriado em {targetdirectory}")
    else:
        os.mkdir(project_path)
        logging.info(f"O diretório {projectname} foi criado em {targetdirectory}")
        print(f"O diretório {projectname} foi criado em {targetdirectory}")

    # Cria subdiretórios
    tables_path = os.path.join(project_path, "tables")
    os.mkdir(tables_path)
    plots_path = os.path.join(project_path, "plots")
    os.mkdir(plots_path)

    return project_path


def read_proteomics_file(rawfilepath, method, control):
    """Lê o arquivo de dados de proteômica e retorna um objeto OmicScope."""
    # Remove espaços e aspas do caminho
    rawfilepath = rawfilepath.strip().strip('"')

    if not os.path.exists(rawfilepath):
        raise FileNotFoundError(
            f"Arquivo não encontrado: {rawfilepath}. Verifique o caminho inserido e tente novamente.")
    elif os.path.getsize(rawfilepath) == 0:
        raise Exception(
            f"O arquivo fornecido está vazio: {rawfilepath}. Verifique o conteúdo do mesmo e tente novamente.")

    try:
        data = omics.OmicScope(rawfilepath, Method=method, ControlGroup=control)
        return data
    except Exception as e:
        logging.error(f"Erro ao ler o arquivo de proteômica: {str(e)}")
        raise Exception(f"Erro ao ler o arquivo de proteômica: {str(e)}")


def plot_data(rawfiledata, analysis_type, titlename, plotsdir, analysis):

    if analysis_type == "conditions_barplot" or analysis_type == "conditions_boxplot":
        proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
        # Remove aspas simples e duplas se o usuário as inserir
        proteins_list = [p.strip().strip("'").strip('"') for p in proteins_string.split(',') if p.strip()]
        palette = input("Insira a paleta desejada (pressione enter para utilizar a paleta de cores padrão): ").lower()
        if len(palette) == 0:
            palette = 'viridis'
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        if analysis_type == "conditions_barplot":
            rawfiledata.bar_protein(*proteins_list,
                                    palette=palette,
                                    dpi=300,
                                    save=os.path.join(plotsdir, titlename))
        else:  # conditions_boxplot
            rawfiledata.boxplot_protein(*proteins_list,
                                 palette=palette,
                                 dpi=300,
                                 save=os.path.join(plotsdir, titlename))

    elif analysis_type == "dynamic_range":
        proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
        proteins_list = [p.strip().strip("'").strip('"') for p in proteins_string.split(',') if p.strip()]
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.DynamicRange(*proteins_list,
                                  dpi=300,
                                  save=os.path.join(plotsdir, titlename))

    elif analysis_type == "ma_plot":
        proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
        proteins_list = [p.strip().strip("'").strip('"') for p in proteins_string.split(',') if p.strip()]
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.MAplot(*proteins_list,
                            dpi=300,
                            save=os.path.join(plotsdir, titlename))

    elif analysis_type == "correlation_heatmap":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.correlation(dpi=300,
                                linewidth=0,
                                save=os.path.join(plotsdir, titlename))

    elif analysis_type == "expression_heatmap":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.heatmap(dpi=300,
                            linewidth=0,
                            save=os.path.join(plotsdir, titlename))

    elif analysis_type == "id_barplot":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.bar_ident(dpi=300,
                           save=os.path.join(plotsdir, titlename))

    elif analysis_type == "pca":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.pca(dpi=300,
                        save=os.path.join(plotsdir, titlename))

    elif analysis_type == "kmeans":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.k_trend(dpi=300,
                           save=os.path.join(plotsdir, titlename))

    elif analysis_type == "normalization_plot":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.normalization_boxplot(dpi=300,
                                       save=os.path.join(plotsdir, titlename))


def perform_ora_enrichment(rawfiledata, db, titlename, analysis, tablesdir, plotsdir):
    print(f"Iniciando análise: {analysis}...")
    logging.info(f"Iniciando análise: {analysis}...")
    data_object = omics.EnrichmentScope(rawfiledata, Analysis='ORA', dbs=[db])
    data_dataframe = data_object.results
    create_xlsx_file(data_dataframe, titlename, tablesdir)
    data_object.dotplot(dpi=300, palette='PuBu', save=os.path.join(plotsdir, f"{titlename}.tiff"))
    print(f"Análise {analysis} concluída.")
    logging.info(f"Análise {analysis} concluída.")


def create_xlsx_file(dataframe, name, tablesdir):
    filename = os.path.join(tablesdir, f"{name}.xlsx")
    dataframe.to_excel(filename)


def filter_deps_dataframe(cutoff, dataframe):
    condition = dataframe['log2(fc)'].abs() > cutoff
    filtered_deps_dataframe = dataframe[condition].copy()
    return filtered_deps_dataframe


def create_json_file(dictionary, name, tablesdir):
    filename = os.path.join(tablesdir, f"{name}.json")
    with open(filename, "w+") as js:
        json.dump(dictionary, js)


# ----------------------------------------------------------------------------------#

if __name__ == "__main__":
    while True:
        try:
            print("Proteoanalyzer 1.0")
            # Gerando diretórios e arquivo de log das análises
            user_name = input("Insira o nome do usuário: ")
            project_name = input("Insira um nome para o projeto a ser analisado: ")
            target_directory = input("Insira o caminho onde os arquivos da análise serão criados (ex: C:/Users/user/Documents): ").strip().strip('"')
            directory = create_dir(target_directory, project_name)
            tables_dir = os.path.join(directory, "tables")
            plots_dir = os.path.join(directory, "plots")

            # Configurando o logger
            log_file_path = os.path.join(directory, f"{project_name}.log")
            logging.basicConfig(level=logging.INFO,
                                format='%(asctime)s - %(levelname)s - %(message)s',
                                filename=log_file_path,
                                filemode='w')

            logging.info(f"Project name: {project_name}")
            logging.info(f"Owned by: {user_name}")
            logging.info(f"Created on: {datetime.now()}")

            print("Diretórios do projeto criados com sucesso.")
            logging.info("Diretórios do projeto criados com sucesso.")

            # Obtendo dados das análises a partir do usuario
            raw_file_path = input(
                "Insira o caminho da planilha com os dados da proteômica (não esqueça de incluir a extensão do arquivo): ").strip()

            # Validação do metodo de proteomica
            valid_methods = ['MaxQuant', 'Progenesis', 'General', 'DIA-NN', 'PatternLab']
            proteomics_method = input("Insira o software de análise dos dados brutos (MaxQuant, Progenesis, General, DIA-NN, PatternLab): ").strip()
            if proteomics_method not in valid_methods:
                print(f"Aviso: Método '{proteomics_method}' pode não ser reconhecido. Métodos válidos: {', '.join(valid_methods)}")
                logging.warning(f"Método de proteômica possivelmente inválido: {proteomics_method}")

            control_group = input("Insira o nome do grupo controle, exatamente como está na planilha: ").strip()

            # Lendo os dados e criando os arquivos com as planilhas
            # Obtendo os dados da proteômica a partir da planilha original
            print("Lendo o arquivo de dados...")
            logging.info("Lendo o arquivo de dados...")
            raw_file_data = read_proteomics_file(raw_file_path, proteomics_method, control_group)
            logging.info("Arquivo de dados lido com sucesso.")

            # Lendo os parÂmetros e salvando-os em um arquivo
            print("Criando tabela de parâmetros...")
            logging.info("Criando tabela de parâmetros...")
            params_dataframe = pd.DataFrame(raw_file_data.Params)
            create_xlsx_file(params_dataframe, "Parâmetros", tables_dir)

            # Criando arquivo json com as condições do estudo
            print("Criando arquivo de condições do estudo...")
            logging.info("Criando arquivo de condições do estudo...")
            conditions_dictionary = {"Condições": f"{raw_file_data.Conditions}",
                                     "Controle": f"{raw_file_data.ControlGroup}"}
            create_json_file(conditions_dictionary, "Condições do estudo", tables_dir)

            # Criando tabela com os dados brutos
            print("Criando tabela de dados brutos...")
            logging.info("Criando tabela de dados brutos...")
            all_data_dataframe = pd.DataFrame(raw_file_data.quant_data)
            create_xlsx_file(all_data_dataframe, "Dados brutos", tables_dir)

            # Criando tabela com as DEPs
            print("Criando tabela de DEPs...")
            logging.info("Criando tabela de DEPs...")
            deps_dataframe = pd.DataFrame(raw_file_data.deps)
            create_xlsx_file(deps_dataframe, "DEPs", tables_dir)

            # Obtendo as DEPs com significância para o estudo
            while True:
                fc_input = input(
                    "Insira o valor de FC que você deseja para calcular o cutoff dos DEPs (1, 1.25, 1.5, 1.75, 2): ")
                if not fc_input:
                    print("Insira um valor de FC para continuar.")
                    logging.warning("Nenhum valor de FC inserido. Por favor, insira um valor para continuar.")
                    continue
                else:
                    break
            try:
                fc = float(fc_input)
                if fc <= 0:
                    raise ValueError("FC deve ser maior que zero.")
                log2fc_cutoff = math.log2(fc)
            except ValueError as e:
                print(f"Erro: Valor de FC inválido. {str(e)}")
                logging.error(f"Valor de FC inválido: {fc_input}")
                continue
            significant_deps_dataframe = filter_deps_dataframe(log2fc_cutoff, deps_dataframe)
            create_xlsx_file(significant_deps_dataframe, "DEPs filtradas", tables_dir)
            logging.info(f"DEPs filtradas com log2(fc) > {log2fc_cutoff:.2f}")

            # Plotando os gráficos
            print("Gerando gráficos...")
            logging.info("Gerando gráficos...")

            # Barplot de identificação
            plot_data(raw_file_data, "id_barplot", f"{project_name}_Barplot_Identificação.tiff", plots_dir,
                      "barplot de identificação")

            # Dynamic range
            plot_data(raw_file_data, "dynamic_range", f"{project_name}_Dynamic_range.tiff", plots_dir,
                      "gráfico de dynamic range")

            # Volcano_plot_selection = input("Deseja plotar o volcano plot? (y/n): ")

            # MA plot
            plot_data(raw_file_data, "ma_plot", f"{project_name}_MA_plot.tiff", plots_dir, "gráfico de MA")

            # Gráfico de normalização
            plot_data(raw_file_data, "normalization_plot", f"{project_name}_Normalização_boxplot.tiff", plots_dir,
                      "gráfico de normalização")

            # Barplot de proteínas entre as condições
            plot_data(raw_file_data, "conditions_barplot", f"{project_name}_Bar_Plot_Proteínas_Entre_Condições.tiff",
                      plots_dir, "gráfico barplot de comparação de proteínas entre condições")

            # Boxplot de proteínas entre as condições
            plot_data(raw_file_data, "conditions_boxplot", f"{project_name}_Boxplot_Proteínas_Entre_Condições.tiff",
                      plots_dir, "gráfico boxplot de comparação de proteínas entre condições")

            # Heatmap de expressão
            plot_data(raw_file_data, "expression_heatmap", f"{project_name}_Heatmap_expressão.tiff", plots_dir,
                      "gráfico heatmap de expressão")

            # Heatmap de correlação
            plot_data(raw_file_data, "correlation_heatmap", f"{project_name}_Heatmap_correlação.tiff", plots_dir,
                      "gráfico heatmap de correlação")

            # PCA
            plot_data(raw_file_data, "pca", f"{project_name}_PCA.tiff", plots_dir, "gráfico de PCA")

            # K-means
            plot_data(raw_file_data, "kmeans", f"{project_name}_Kmeans.tiff", plots_dir, "gráfico de K-means")

            print("Plotagem dos dados concluída.")
            logging.info("Plotagem dos dados concluída.")

            # Enriquecimento dos dados
            print("Iniciando enriquecimento dos dados...")
            logging.info("Iniciando enriquecimento dos dados...")

            # KEGG
            perform_ora_enrichment(raw_file_data, 'KEGG_2021_Human', f"{project_name}_KEGG_pathways", "vias KEGG",
                                   tables_dir, plots_dir)

            # GO: processo biológico
            perform_ora_enrichment(raw_file_data, 'GO_Biological_Process_2025', f"{project_name}_GO_BP","processo biológico", tables_dir, plots_dir)

            # CC: componente celular
            perform_ora_enrichment(raw_file_data, 'GO_Cellular_Component_2025', f"{project_name}_GO_CC", "componente celular", tables_dir, plots_dir)

            # MF: função molecular
            perform_ora_enrichment(raw_file_data, 'GO_Molecular_Function_2025', f"{project_name}_GO_MF", "função molecular", tables_dir, plots_dir)

            # Reactome
            perform_ora_enrichment(raw_file_data, 'Reactome_Pathways_2024', f"{project_name}_Reactome", "vias Reactome", tables_dir, plots_dir)

            # OMIM
            perform_ora_enrichment(raw_file_data, 'OMIM_Expanded', f"{project_name}_OMIM", "vias OMIM", tables_dir, plots_dir)

            # DisGeNET
            perform_ora_enrichment(raw_file_data, 'DisGeNET', f"{project_name}_DisGeNET", "vias DisGeNET", tables_dir, plots_dir)
            print("Enriquecimento dos dados concluído.")
            logging.info("Enriquecimento dos dados concluído.")
            print("Análise dos dados concluída com sucesso!")
            logging.info("Análise dos dados concluída com sucesso!")

            while True:
                repeat = input("Deseja realizar uma nova análise? (y/n): ").lower().strip()
                if repeat == "y":
                    break  # Sai do loop interno e continua o loop principal
                elif repeat == "n":
                    print("Encerrando o programa...")
                    logging.info("Programa encerrado pelo usuário.")
                    exit(0)  # Encerra o programa
                else:
                    print(f"Valor inválido inserido: '{repeat}'. Por favor, digite 'y' para sim ou 'n' para não.")
                    logging.warning(f"Valor inválido inserido: {repeat}")

        except ValueError as error:
            print(f"O valor inserido é invalido: {error}")
            logging.error(f"O valor inserido é invalido: {error}", exc_info=True)
        except Exception as error:
            print(f"Erro inesperado: {error}")
            logging.error(f"Erro inesperado: {error}", exc_info=True)









