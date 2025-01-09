import os
import PyPDF2
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ProcessPoolExecutor, as_completed

def plot_protein_boxplot(data, output_dir):
    protein_id, group_data = data
    try:
        file_path = os.path.join(output_dir, f'{protein_id}.pdf')
        with PdfPages(file_path) as pdf:
            plt.figure(figsize=(8, 6))
            group_data.boxplot(column=protein_id, by='group', grid=False)
            plt.title(f'Boxplot of Protein Expression for {protein_id}')
            plt.xlabel('')
            plt.ylabel('Expression Level')
            pdf.savefig()
            plt.close()
        return file_path
    except Exception as e:
        print(f"Error processing {protein_id}: {e}")
        return None


def merge_pdfs(files, output_dir):
    pdf_writer = PyPDF2.PdfWriter()
    for file in files:
        if file is not None and os.path.exists(file):
            pdf_reader = PyPDF2.PdfReader(file)
            for page in range(len(pdf_reader.pages)):
                pdf_writer.add_page(pdf_reader.pages[page])
            os.remove(file)
    output_pdf_path = os.path.join(output_dir, 'qc_boxplots.pdf')
    with open(output_pdf_path, 'wb') as out:
        pdf_writer.write(out)
    print(f'Saved combined plots to {output_pdf_path}')


def parallel_plot(dataframe, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    protein_ids = dataframe.columns.drop('group')
    files = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        future_to_protein = {executor.submit(plot_protein_boxplot, (pid, dataframe[['group', pid]]), output_dir): pid for pid in protein_ids}
        for future in as_completed(future_to_protein):
            result = future.result()
            if result:
                files.append(result)

    merge_pdfs(files, output_dir)



