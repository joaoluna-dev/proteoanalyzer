# Created by João Gabriel Barbosa de Luna - iLIKA/UFPE
import json
import omicscope as omics
import pandas as pd
import seaborn as sns
import numpy as np
import os
from datetime import datetime
import math
import logging
import shutil
from matplotlib import pyplot as plt


# ---------------------------------------------------------------------------------#
def create_dir(targetdirectory, projectname):
    """Cria a estrutura de diretórios do projeto."""
    if not os.path.exists(targetdirectory):
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
                    if os.path.exists(project_path):
                        print(f"O diretório {projectname} também existe no caminho selecionado.")
                        continue
                    else:
                        print(f"O diretório {projectname} foi criado em {another_directory}")
                        break
                elif another_dir == "n":
                    print("Operação cancelada pelo usuário.")
                    exit(0)
                else:
                    print(f"Opção inválida inserida: {another_dir}")
                    continue
        elif response == 'y':
            shutil.rmtree(project_path)
            os.mkdir(project_path)
            print(f"O diretório {projectname} foi recriado em {targetdirectory}")
    else:
        os.mkdir(project_path)
        print(f"O diretório {projectname} foi criado em {targetdirectory}")

    # Cria subdiretórios e arquivo de log
    tables_path = os.path.join(project_path, "tables")
    os.mkdir(tables_path)
    plots_path = os.path.join(project_path, "plots")
    os.mkdir(plots_path)

    # Criando sistema de log
    log_file_path = os.path.join(project_path, "app.log")
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filename=log_file_path,
                        filemode='w')

    logging.info(f"Project name: {project_name}")
    logging.info(f"Owned by: {user_name}")
    logging.info(f"Created on: {datetime.now()}")
    logging.info(f"O diretório {project_path} foi criado.")

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
                                    save=os.path.join(plotsdir, f"{titlename}_"))
        else:  # conditions_boxplot
            rawfiledata.boxplot_protein(*proteins_list,
                                 palette=palette,
                                 dpi=300,
                                 save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "dynamic_range":
        proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
        proteins_list = [p.strip().strip("'").strip('"') for p in proteins_string.split(',') if p.strip()]
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.DynamicRange(*proteins_list,
                                  dpi=300,
                                  save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "ma_plot":
        proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
        proteins_list = [p.strip().strip("'").strip('"') for p in proteins_string.split(',') if p.strip()]
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.MAplot(*proteins_list,
                            dpi=300,
                            save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "correlation_heatmap":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.correlation(dpi=300,
                                linewidth=0,
                                save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "expression_heatmap":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.heatmap(dpi=300,
                            linewidth=0,
                            save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "id_barplot":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.bar_ident(dpi=300,
                           save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "pca":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.pca(dpi=300,
                        save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "kmeans":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.k_trend(dpi=300,
                           save=os.path.join(plotsdir, f"{titlename}_"))

    elif analysis_type == "normalization_plot":
        print(f"Gerando {analysis}...")
        logging.info(f"Gerando {analysis}...")
        rawfiledata.normalization_boxplot(dpi=300,
                                       save=os.path.join(plotsdir, f"{titlename}_"))

def create_volcano_plot(df, cutoff, pvalue_cutoff, plotsdir, titlename):
    proteins_string = input("Insira as proteínas que deseja apontar no gráfico, separadas por vírgula. Certifique-se de que as proteínas indicadas estão na tabela, e que estão escritas exatamente como estão na tabela: ")
    proteins_list = [p.strip() for p in proteins_string.split(',') if p.strip()]
    full_path = os.path.join(plotsdir, f"{titlename}.svg")
    df['regulation'] = 'Not Significant'
    condition_up = (df['log2(fc)'] >= cutoff) & (df['pAdjusted'] < pvalue_cutoff)
    condition_down = (df['log2(fc)'] <= -cutoff) & (df['pAdjusted'] < pvalue_cutoff)
    df.loc[condition_up, 'regulation'] = 'Up-regulated'
    df.loc[condition_down, 'regulation'] = 'Down-regulated'

    # --- Plot Styling ---
    plt.style.use('seaborn-v0_8-white')
    plt.figure(figsize=(8, 8))

    palette = {
        'Up-regulated': '#ff4d4d',
        'Down-regulated': '#4d4dff',
        'Not Significant': 'darkgrey'
    }

    sns.scatterplot(
        data=df, x='log2(fc)', y='-log10(pAdjusted)', hue='regulation',
        palette=palette, s=40, alpha=0.7, edgecolor='none'
    )

    # --- Axis and threshold lines ---
    plt.axvline(x=0, color='black', linestyle='-', linewidth=0.75)
    plt.axhline(y=-np.log10(pvalue_cutoff), color='dimgrey', linestyle='--', linewidth=1)
    plt.axvline(x=cutoff, color='dimgrey', linestyle='--', linewidth=1)
    plt.axvline(x=-cutoff, color='dimgrey', linestyle='--', linewidth=1)

    # --- Protein Labeling ---
    df_labels = df[df['gene_name'].isin(proteins_list)]
    for index, row in df_labels.iterrows():
        plt.text(row['log2(fc)'], row['-log10(pAdjusted)'], row['gene_name'],
                 ha='center', va='bottom', fontsize=9, fontweight='bold', alpha=0.9)

    # --- Final Plot Customization ---
    plt.title('Volcano Plot', fontsize=16, fontweight='bold')
    plt.ylabel(r'$-log_{10}(Adjusted\ p-value)$', fontsize=12)
    plt.xlabel(r'$log_2(Fold\ Change)$', fontsize=12)
    plt.xlim(-3, 3)

    handles, labels = plt.gca().get_legend_handles_labels()
    order = ['Up-regulated', 'Down-regulated', 'Not Significant']
    plt.legend([handles[labels.index(key)] for key in order], order, title='Regulation', frameon=False)

    try:
        print("Gerando volcano plot...")
        plt.savefig(full_path, format='svg', dpi=300, bbox_inches='tight')
        print(f"\nVolcano plot salvo com sucesso em: '{full_path}'")
        logging.info(f"Volcano plot salvo com sucesso em: '{full_path}'")
    except Exception as e:
        print(f"Um erro ocorreu durante a geração do volcano plot: {e}")
        logging.error(f"Um erro ocorreu durante a geração do volcano plot: {e}")


def perform_ora_enrichment(rawfiledata, db, titlename, analysis, tablesdir, plotsdir):
    print(f"Iniciando análise de: {analysis}...")
    logging.info(f"Iniciando análise: {analysis}...")
    data_object = omics.EnrichmentScope(rawfiledata, Analysis='ORA', dbs=[db])
    data_dataframe = data_object.results
    create_xlsx_file(data_dataframe, titlename, tablesdir)
    data_object.dotplot(dpi=300, palette='PuBu', save=os.path.join(plotsdir, f"{titlename}_"))
    print(f"Análise de {analysis} concluída.")
    logging.info(f"Análise de {analysis} concluída.")


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


            print("Diretórios do projeto criados com sucesso.")
            logging.info("Diretórios do projeto criados com sucesso.")

            # Obtendo dados das análises a partir do usuario
            raw_file_path = input("Insira o caminho da planilha com os dados da proteômica (não esqueça de incluir a extensão do arquivo): ").strip()

            # Validação do metodo de proteomica
            valid_methods = ['MaxQuant', 'Progenesis', 'General', 'DIA-NN', 'PatternLab']
            while True:
                proteomics_method = input("Insira o software de análise dos dados brutos (MaxQuant, Progenesis, General, DIA-NN, PatternLab): ").strip()
                if proteomics_method not in valid_methods:
                    print(f"Aviso: Método '{proteomics_method}' pode não ser reconhecido. Métodos válidos: {', '.join(valid_methods)}. Tente novamente.")
                    logging.warning(f"Método de proteômica possivelmente inválido: {proteomics_method}")
                    continue
                else:
                    break

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
            print("Iniciando plotagem de barplot de identificação das condições...")
            plot_data(raw_file_data, "id_barplot", project_name, plots_dir,
                      "barplot de identificação")

            # Dynamic range
            print("Iniciando plotagem de Dynamic Range...")
            plot_data(raw_file_data, "dynamic_range", project_name, plots_dir,
                      "gráfico de dynamic range")

            # Volcano plot
            print("Iniciando plotagem de Volcano Plot...")
            create_volcano_plot(deps_dataframe, log2fc_cutoff, 0.05, plots_dir, f"{project_name}_Volcano_plot")

            # MA plot
            print("Iniciando plotagem de MA plot...")
            plot_data(raw_file_data, "ma_plot", project_name, plots_dir, "gráfico de MA")

            # Gráfico de normalização
            print("Iniciando plotagem de gráfico de normalização...")
            plot_data(raw_file_data, "normalization_plot", project_name, plots_dir,
                      "gráfico de normalização")

            # Barplot de proteínas entre as condições
            print("Iniciando plotagem de barplot comparando proteínas entre condições...")
            plot_data(raw_file_data, "conditions_barplot", project_name,
                      plots_dir, "gráfico barplot de comparação de proteínas entre condições")

            # Boxplot de proteínas entre as condições
            print("Iniciando plotagem de boxplot comparando proteínas entre condições...")
            plot_data(raw_file_data, "conditions_boxplot", project_name,
                      plots_dir, "gráfico boxplot de comparação de proteínas entre condições")

            # Heatmap de expressão
            print("Iniciando plotagem de heatmap de expressão...")
            plot_data(raw_file_data, "expression_heatmap", project_name, plots_dir,
                      "gráfico heatmap de expressão")

            # Heatmap de correlação
            print("Iniciando plotagem de heatmap de correlação...")
            plot_data(raw_file_data, "correlation_heatmap", project_name, plots_dir,
                      "gráfico heatmap de correlação")

            # PCA
            print("Iniciando plotagem de gráfico de PCA...")
            plot_data(raw_file_data, "pca", project_name, plots_dir, "gráfico de PCA")

            # K-means
            print("Iniciando plotagem de gráfico de PCA...")
            plot_data(raw_file_data, "kmeans", project_name, plots_dir, "gráfico de K-means")

            print("Plotagem dos dados concluída.")
            logging.info("Plotagem dos dados concluída.")

            # Enriquecimento dos dados
            print("Iniciando enriquecimento dos dados...")
            logging.info("Iniciando enriquecimento dos dados...")

            # KEGG
            perform_ora_enrichment(raw_file_data, 'KEGG_2021_Human', project_name, "vias KEGG",
                                   tables_dir, plots_dir)

            # GO: processo biológico
            perform_ora_enrichment(raw_file_data, 'GO_Biological_Process_2025', project_name,"processo biológico", tables_dir, plots_dir)

            # CC: componente celular
            perform_ora_enrichment(raw_file_data, 'GO_Cellular_Component_2025', project_name, "componente celular", tables_dir, plots_dir)

            # MF: função molecular
            perform_ora_enrichment(raw_file_data, 'GO_Molecular_Function_2025', project_name, "função molecular", tables_dir, plots_dir)

            # Reactome
            perform_ora_enrichment(raw_file_data, 'Reactome_Pathways_2024', project_name, "vias Reactome", tables_dir, plots_dir)

            # OMIM
            perform_ora_enrichment(raw_file_data, 'OMIM_Expanded', project_name, "vias OMIM", tables_dir, plots_dir)

            # DisGeNET
            perform_ora_enrichment(raw_file_data, 'DisGeNET', project_name, "vias DisGeNET", tables_dir, plots_dir)
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









