# Mol_Similarity_GCMSD

This is an application created in python libraries like RDKIT, plotly, streamlet, where is possible to analyze the data of the aroma of coffee samples obtained  by Gas chromatography–mass spectrometry. The data was obtained from a thesis work called "characterization of the aroma of ground coffee from Puerto Rico using the technique of solid phase microextraction (spme) and gas chromatography coupled to mass spectrometry (gc/ms)" (link = https://avdiaz.files.wordpress.com/2010/09/rojasmonroy.pdf).
In this streamlit application it is possible to load the data obtained from the chemical analysis of two coffee samples and calculate the molecular similarity between them, molecule by molecule with a 3d representation or by batch of molecules per retention time range.
The plotly plot shows the similarity by three regions of retention time 0 to 20 min, 20 to 40 min and 40 to 60 min. The first segment has the most volatile molecules, in the second segment there are fewer molecule volatility, and the third part, the molecules have the fewest vapor pressure. With this information, we can compare the aroma performance of coffee samples over the time.




![similarity_GCMSD_2](https://user-images.githubusercontent.com/59380458/229645710-ab074425-c0cf-4136-a95c-2e00a457f2a3.gif)
