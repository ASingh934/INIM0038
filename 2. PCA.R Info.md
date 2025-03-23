**PCA.R**

PCA is performed on _**gene_data**_. For each PC, the Standard Deviation and proportion of total variance is displayed. For PC1, the top 50 most weighted genes are displayed to show which genes are most effective at clustering the patients. Then, 2 PCA plots are drawn.

**PCA plot 1**: All patients are plotted here, grouped according to whether they are BioAID or controls. We can see from visual inspection that the controls cluster in the Right.

![image](https://github.com/user-attachments/assets/314b2404-061b-4ead-b66f-267eedd52f9c)

**PCA plot 2**: For the Oxford & UCL patients whose clinical labels are available, they are replotted on the same PCA plot, and grouped by whether or not they have an ARI (and if so, then the pathogen which is causing it). From this, we can see no clear clustering of any group in particular, and we conclude that PCA plots are not an effective method of endotyping ARI patients.

![image](https://github.com/user-attachments/assets/11d5d1af-ee41-4d3b-a153-d01b4f8f8de5)

A scree plot is shown in the code for completeness.
