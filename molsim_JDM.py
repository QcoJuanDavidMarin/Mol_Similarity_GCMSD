import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import glob

from rdkit.Chem import AllChem as Chem
import rdkit
from rdkit import DataStructs
from rdkit.Chem import Draw, PandasTools, AllChem

import plotly.express as px
import numpy as np
import seaborn as sns
import warnings
import streamlit as st
import matplotlib.pyplot as plt
from io import StringIO
import py3Dmol

import warnings
warnings.filterwarnings('ignore')


st.set_page_config(page_title='Molecular_similarity_CGMSD_in_Coffee_JDM',
                    page_icon='👁',layout="wide")

def main():

    st.title('A streamlit app to calculate :red[TANIMOTO SIMILARITY] of molecules identified by _Gas chromatography–mass spectrometry (GC-MS)_ in :blue[_coffee samples_]')

    st.success(':red[Description]: This is an application created in python libraries like RDKIT, plotly, streamlet, where is possible to analyze the data of the aroma of coffee samples obtained  by Gas chromatography–mass spectrometry. The data was obtained from a thesis work called "characterization of the aroma of ground coffee from Puerto Rico using the technique of solid phase microextraction (spme) and gas chromatography coupled to mass spectrometry (gc/ms)" (link = https://avdiaz.files.wordpress.com/2010/09/rojasmonroy.pdf). In this streamlit application it is possible to load the data obtained from the chemical analysis of two coffee samples and calculate the molecular similarity between them, molecule by molecule with a 3d representation or by batch of molecules per retention time range. The plotly plot shows the similarity by three regions of retention time 0 to 20 min, 20 to 40 min and 40 to 60 min. The first segment has the most volatile molecules, in the second segment there are fewer molecule volatility, and the third part, the molecules have the fewest vapor pressure. With this information, we can compare the aroma performance of coffee samples over the time.')
    



    st.subheader(':blue[Streamlit app Developed by:] **Juan David Marín**')
    st.markdown("LinkInd: [link](https://www.linkedin.com/in/qco-juan-david-marin/)")
    st.markdown("Github: [link](https://github.com/QcoJuanDavidMarin)")

    # Create a file uploader widget that accepts multiple CSV files
    uploaded_files_CGMDS_1 = st.sidebar.file_uploader('Upload GCMSD file_1', type='csv')
    if uploaded_files_CGMDS_1 is not None:
        file_CGMDS_1= pd.read_csv(uploaded_files_CGMDS_1, sep=";", encoding='cp1252')
        file_CGMDS_1['place'] = uploaded_files_CGMDS_1.name

    uploaded_files_CGMDS_2 = st.sidebar.file_uploader('Upload GCMSD file_2', type='csv')
    if uploaded_files_CGMDS_2 is not None:
        file_CGMDS_2= pd.read_csv(uploaded_files_CGMDS_2, sep=";", encoding='cp1252')
        file_CGMDS_2['place'] = uploaded_files_CGMDS_2.name

    df_coffe_cgmsd = pd.concat([file_CGMDS_1,file_CGMDS_2],axis = 0, ignore_index= True)
    #st.dataframe(df_coffe_cgmsd)
    coffee_city_1 = df_coffe_cgmsd[df_coffe_cgmsd['place'] == uploaded_files_CGMDS_1.name].reset_index(drop=True)
    coffee_city_2 = df_coffe_cgmsd[df_coffe_cgmsd['place'] == uploaded_files_CGMDS_2.name].reset_index(drop=True)

    # Calculing data frame similarity 
    slider_tr = st.sidebar.slider('Select time retention (min)', 0, 60, (5, 10))

    #st.dataframe(df_coffe_cgmsd)


    ##_________________________________________
    ## Molecles in 3D and calculing similarity
    col1, col2 = st.columns(2)
    with col1:
        number_1 = st.number_input('Choose molecule CGMSD_1',0,60)
        def show(smi, style='stick'):
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            #AllChem.MMFFOptimizeMolecule(mol, maxIters = 200)
            mblock = Chem.MolToMolBlock(mol)

            view = py3Dmol.view(width=400, height=400)
            view.addModel(mblock, 'mol')
            view.setStyle({style:{}})
            view.addSurface('VDW', {'opacity':0.5, 'colorscheme':{'gradient':'rwb'}})
            view.zoomTo()
            view.show()
            view.render()
            t =view.js()
            f = open('viz.html', 'w')
            f.write(t.startjs)
            f.write(t.endjs)
            f.close()
        compound_smiles_1 = coffee_city_1['smiles'][number_1]
        m = Chem.MolFromSmiles(compound_smiles_1)
        Draw.MolToFile(m,'mol.png')
        show(compound_smiles_1)
        HtmlFile = open("viz.html", 'r', encoding='utf-8')
        source_code1 = HtmlFile.read()
        HtmlFile = open("viz.html", 'r', encoding='utf-8')
        source_code1 = HtmlFile.read()
        components.html(source_code1, height = 400,width=400)
        st.write(coffee_city_1['Compound'][number_1])

    with col2:
        number_2 = st.number_input('Choose molecule CGMSD_2',0,60)
        def show(smi, style='stick'):
            mol2 = Chem.MolFromSmiles(smi)
            mol2 = Chem.AddHs(mol2)
            AllChem.EmbedMolecule(mol2)
            #AllChem.MMFFOptimizeMolecule(mol, maxIters = 200)
            mblock = Chem.MolToMolBlock(mol2)
            view = py3Dmol.view(width=400, height=400)
            view.addModel(mblock, 'mol')
            view.setStyle({style:{}})
            view.addSurface('VDW', {'opacity':0.5, 'colorscheme':{'gradient':'rwb'}})
            view.zoomTo()
            view.show()
            view.render()
            t =view.js()
            f = open('viz.html', 'w')
            f.write(t.startjs)
            f.write(t.endjs)
            f.close()
        compound_smiles_2 = coffee_city_2['smiles'][number_2]
        m = Chem.MolFromSmiles(compound_smiles_2)
        Draw.MolToFile(m,'mol.png')
        show(compound_smiles_2)
        HtmlFile = open("viz.html", 'r', encoding='utf-8')
        source_code2 = HtmlFile.read()
        HtmlFile = open("viz.html", 'r', encoding='utf-8')
        source_code2 = HtmlFile.read()
        components.html(source_code2, height = 400,width=400)
        st.write(coffee_city_2['Compound'][number_2])
    ## Similarity molecules 3D
    css='''
    [data-testid="metric-container"] {
        width: fit-content;
        margin: auto;
    }
    [data-testid="metric-container"] > div {
        width: fit-content;
        margin: auto;
    }
    [data-testid="metric-container"] label {
        width: fit-content;
        margin: auto;
    }
    '''
    st.markdown(f'<style>{css}</style>',unsafe_allow_html=True)
    mfp1 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(compound_smiles_1), radius=2,nBits = 1024)
    mfp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(compound_smiles_2), radius=2,nBits = 1024)
    TaniSimil = round(DataStructs.TanimotoSimilarity(mfp1, mfp2),3)
    st.metric(label = 'Tanimoto Similarity',value=TaniSimil, delta = f"CGMSD_1: {coffee_city_1['Compound'][number_1]} Vs CGMSD_2: {coffee_city_2['Compound'][number_2]}")


    ## images of molecules in 2D 
    cola, colb = st.columns(2, gap = "large")
    with cola:
        ## Rendering molecules 2D data 1
        def show_mol_1(df1,df2,rt_star, rt_end):
            df11 = df1.copy(deep=True)
            df22 = df2.copy(deep=True)
            PandasTools.AddMoleculeColumnToFrame(df11,'smiles','ROMol1')
            PandasTools.AddMoleculeColumnToFrame(df22,'smiles','ROMol2')
            # Concatenate both dataframes
            df_both1 = pd.concat([df11, df22], axis = 1).dropna().reset_index(drop=True)
            ## Add a number sequence to identify each column with the same name
            coluns = [str(i) + '_' + col for i, col in enumerate(df_both1.columns)]
            df_both1.columns = coluns
            df_both1 = df_both1[(df_both1['19_retention_time'] >= rt_star) & (df_both1['19_retention_time'] <= rt_end)]

            rMols1 = [Chem.MolFromSmiles(r) for r in df_both1['6_smiles']]
            rMols1 = Draw.MolsToGridImage(rMols1,legends=[label for label in df_both1['0_Compound']],
                        subImgSize=(150,150), useSVG=False, molsPerRow= 7)
            rMols1 = st.image(rMols1, caption='Molecules DataFrame 1')
            return(rMols1)
        show_mol_1(coffee_city_1,coffee_city_2,slider_tr[0],slider_tr[1])
    with colb:
            ## Rendering molecules 2D data 2
        def show_mol_2(df1,df2,rt_star, rt_end):
            df11 = df1.copy(deep=True)
            df22 = df2.copy(deep=True)
            PandasTools.AddMoleculeColumnToFrame(df11,'smiles','ROMol1')
            PandasTools.AddMoleculeColumnToFrame(df22,'smiles','ROMol2')
            # Concatenate both dataframes
            df_both1 = pd.concat([df11, df22], axis = 1).dropna().reset_index(drop=True)
            ## Add a number sequence to identify each column with the same name
            coluns = [str(i) + '_' + col for i, col in enumerate(df_both1.columns)]
            df_both1.columns = coluns
            df_both1 = df_both1[(df_both1['19_retention_time'] >= rt_star) & (df_both1['19_retention_time'] <= rt_end)]

            rMols2 = [Chem.MolFromSmiles(r) for r in df_both1['17_smiles']]
            rMols2 = Draw.MolsToGridImage(rMols2,legends=[label for label in df_both1['11_Compound']],
                        subImgSize=(150,150), useSVG=False, molsPerRow= 7)
            rMols2 = st.image(rMols2, caption='Molecules DataFrame 2')
            return(rMols2)
        show_mol_2(coffee_city_1,coffee_city_2,slider_tr[0],slider_tr[1])
    ##_________________________________________
    ## Function to calculate heatmap
    def simil_matrix(df1,df2, rt_star, rt_end):
        # Display all dataset
        pd.options.display.max_columns = None
        pd.options.display.max_rows = None
        df11 = df1.copy(deep=True)
        df22 = df2.copy(deep=True)
        PandasTools.AddMoleculeColumnToFrame(df11,'smiles','ROMol1')
        PandasTools.AddMoleculeColumnToFrame(df22,'smiles','ROMol2')
        # create de morganfingerprint as vector
        df11['mfp1'] = df11['smiles'].apply(lambda x:Chem.MolFromSmiles(x)).apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2,nBits = 1024))
        df22['mfp2'] = df22['smiles'].apply(lambda x:Chem.MolFromSmiles(x)).apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2,nBits = 1024))
        # Concatenate both dataframes
        df_both1 = pd.concat([df11, df22], axis = 1).dropna().reset_index(drop=True)
        ## Add a number sequence to identify each column with the same name
        coluns = [str(i) + '_' + col for i, col in enumerate(df_both1.columns)]
        df_both1.columns = coluns
        
        df_both = df_both1[(df_both1['20_retention_time'] >= rt_star) & (df_both1['20_retention_time'] <= rt_end)]

        for i in df_both.index:
            fp1 = df_both.loc[i,'11_mfp1']
            colname = df_both.loc[i,'0_Compound']
            sim_matrix = []
            for mol in df_both['22_ROMol2']:
                fp = Chem.GetMorganFingerprintAsBitVect(mol, radius = 2, nBits = 1024)
                sim = round(DataStructs.TanimotoSimilarity(fp1, fp),2)
                sim_matrix.append(sim)
            df_both[colname] = sim_matrix

        #images = list(df_both['10_ROMol1'])
        headers = list(df_both.iloc[:,24:].columns)
        #headers = images
        x = pd.DataFrame(df_both.iloc[:,24:]) 
        x.columns = headers 

        x2 = df_both[['20_retention_time','12_Compound']]
        #x2.columns = [x2.columns,['12_Compound', '20_retention_time', '22_ROMol2']]
        df_final = pd.concat([x2, x], axis=1).reset_index(drop=True)
        cm = sns.light_palette("orange", as_cmap=True)
        #df_final = df_final.style.background_gradient(cmap=cm)
        return(st.write(df_final.style.background_gradient(cmap=cm)))
    simil_matrix(coffee_city_1,coffee_city_2, slider_tr[0],slider_tr[1])
    ##_________________________________________
    ##_________________________________________
    #plotting
    def simil_plot(df1,df2):
            # Display all dataset
        df11 = df1.copy(deep=True)
        df22 = df2.copy(deep=True)
        PandasTools.AddMoleculeColumnToFrame(df11,'smiles','ROMol1')
        PandasTools.AddMoleculeColumnToFrame(df22,'smiles','ROMol2')
            # create de morganfingerprint as vector
        df11['mfp1'] = df11['smiles'].apply(lambda x:Chem.MolFromSmiles(x)).apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2,nBits = 1024))
        df22['mfp2'] = df22['smiles'].apply(lambda x:Chem.MolFromSmiles(x)).apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2,nBits = 1024))
            # Concatenate both dataframes
        df_both1 = pd.concat([df11, df22], axis = 1).dropna().reset_index(drop=True)
            ## Add a number sequence to identify each column with the same name
        coluns = [str(i) + '_' + col for i, col in enumerate(df_both1.columns)]
        df_both1.columns = coluns
        # plotting
        df_both1['TanimotoSimilarity'] = list(map(lambda x,y: DataStructs.TanimotoSimilarity(x,y), df_both1['11_mfp1'], df_both1['23_mfp2']))
        conditionlist = [(df_both1['20_retention_time'] < 20) ,(df_both1['20_retention_time'] >= 20) & (df_both1['8_retention_time'] < 40),(df_both1['8_retention_time'] >= 40)]
        choicelist = ['0-20 min', '20-40 min', '40-60 min']
        df_both1['RT_Range'] = np.select(conditionlist, choicelist, default='Not Specified')
        fig = px.scatter(df_both1, x = '20_retention_time', y = 'TanimotoSimilarity', color='RT_Range',
                                marginal_y = 'box', marginal_x='histogram',template = "simple_white",trendline="ols",
                                title='Comparation between Coffe', hover_data=['0_Compound','12_Compound'])
        fig.update_layout(xaxis_title='RT min', yaxis_title='Tanimoto Similarity')
        # trace_contour = px.density_contour(df_both1, x = '8_retention_time', y = 'TanimotoSimilarity', color='RT_Range',).update_layout(showlegend=False).select_traces()
        # fig.add_traces(list(trace_contour))
        fig.update_layout(legend_title_text="Cromatography time")
        fig.update_xaxes(tickvals=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60])
        fig.update_xaxes(rangeslider_visible=False)
        fig.update_layout(autosize=False,width=1200,height=600)
        return st.plotly_chart(fig,use_container_width=False)
    simil_plot(coffee_city_1,coffee_city_2)


    ## plot of vapor pressure
    with st.expander('**Show vapor pressure**'):
        firstcol, secondcol = st.columns(2,gap = "large")
        with firstcol:
            def simil_matrix(df1):
                # Display all dataset
                df11 = df1.copy(deep=True)
                conditionlist = [(df11['retention_time'] < 20) ,(df11['retention_time'] >= 20) & (df11['retention_time'] < 40),(df11['retention_time'] >= 40)]
                choicelist = ['0-20 min', '20-40 min', '40-60 min']
                df11['RT_Range'] = np.select(conditionlist, choicelist, default='Not Specified')
                fig_1 = px.scatter(df11, x = 'retention_time', y = 'Vapor_pressure', color='RT_Range',
                                        template = "simple_white",log_y=True,
                                        title='Vapor pressure data GCMSD coffee 1')
                fig_1.update_layout(xaxis_title='RT min', yaxis_title='Vapor pressure')
                fig_1.update_layout(autosize=False,width=500,height=400)
                return st.plotly_chart(fig_1,use_container_width=False)
            simil_matrix(coffee_city_1)
        with secondcol:
            def simil_matrix(df2):
                # Display all dataset
                df22 = df2.copy(deep=True)
                conditionlist = [(df22['retention_time'] < 20) ,(df22['retention_time'] >= 20) & (df22['retention_time'] < 40),(df22['retention_time'] >= 40)]
                choicelist = ['0-20 min', '20-40 min', '40-60 min']
                df22['RT_Range'] = np.select(conditionlist, choicelist, default='Not Specified')
                fig_2 = px.scatter(df22, x = 'retention_time', y = 'Vapor_pressure', color='RT_Range',
                                        template = "simple_white",log_y=True,
                                        title='Vapor pressure data GCMSD coffee 2')
                fig_2.update_layout(autosize=False,width=500,height=400)
                fig_2.update_layout(xaxis_title='RT min', yaxis_title='Vapor pressure')
                return st.plotly_chart(fig_2,use_container_width=False)
            simil_matrix(coffee_city_2)

    ##_________________________________________




if __name__ == '__main__':
	main()
