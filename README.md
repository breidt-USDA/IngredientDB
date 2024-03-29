# IngredientDB
Buffer modeling software for formulation and safety of acid/acidified foods

IngredientDB 
A graphical user interface application for estimating the pH of acid and acidified foods based on food ingredient buffer models. 
Software downloads available: https://www.ars.usda.gov/southeast-area/raleigh-nc/fsmqhru/
Fred Breidt, PhD, USDA/ARS Research Microbiologist 
Food Science Market Quality and Handling Research Unit, 322 Schaub Hall, Box 7624, NC State University, Raleigh, NC. 

Overview: Predictions of the final pH of mixtures of food ingredients with food-grade acids such as vinegar have been difficult because each food ingredient may influence pH differently. Quantitatively measuring how pH will change in a food may be helpful for predicting the safety and microbial stability of the food. To aid in this process IngredientDB and a database with information on how foods influence pH have been developed. The program uses database files generated from titration data using previously published USDA/ARS software BufferCapacity3 (reference 1 below). For further details of buffer modeling see references 2-4 below.  

1.	Breidt F. 2023. BufferCapacity3 an interactive GUI program for modelling food ingredient buffering and pH. SoftwareX. 22:101351. https://doi.org/10.1016/j.softx.2023.101351)
2.	Longtin M, Price RE, Mishra R, Breidt F. 2020. Modeling the buffer capacity of ingredients in salad dressing products. J Food Sci 85(4):910-917. https://doi.org/10.1111/1750-3841.15018
3.	Price RE, Longtin M, Conley Payton S, Osborne JA, Johanningsmeier SD, Bitzer D, Breidt F. 2020. Modeling buffer capacity and pH in acid and acidified foods. J Food Sci 85(4):918-925. https://doi.org/10.1111/1750-3841.15091
4.	Breidt F, Skinner CR. 2022. Buffer models for pH and acid changes occurring in cucumber juice fermented with Lactiplantibacillus pentosus and Leuconostoc mesenteroides. J Food Prot 85(9):1273-1281. https://doi.org/10.4315/JFP-22-068

One or more database files generated by BufferCapacity3 may be opened in IngredientDB to create in silico ingredient mixtures. The estimated pH of the ingredient mixture and the percent buffering attributed to each ingredient is then reported in a table. Selected acids can then be added to the mixture, updating the pH and buffering information. Formulations can then be saved as a new database file for future use. The data on pH and buffering percentages can be useful for product development and can help estimate pH stability of the formulation by changing ingredient concentrations. 

Steps for the analysis of ingredient formulations:
1)	Open database file(s): The file | Open Database File(s) menu item allows the user to select one or more â€˜.matâ€™ database files as generated by the BufferCapacity3 software from titration data (see details above)
2)	Adjust ingredient concentrations as needed: Double clicking on the ingredient concentration will allow the user to change the concentration (in percent) of individual ingredients in a formulation. This will automatically update the pH and other data values as well as the Formulation Buffering Percentages table.
  -	Note that the NaCl% value is based on the ionic strength value (from NaCl percent) from the initial titration. This value is automatically set.  
  -	Ingredients may be deleted from the Ingredient Table by double clicking on the ingredient name in the table and responding to the confirmation message box affirmatively.  
3)	Add acid(s) from the Acid table: Acids can be added to the formulation by changing the concentration value either as millimolar (mM) or percent concentration. If the salt of an acid is used (i.e. sodium acetate) clicking on the (Salt) checkbox for a given acid will update the pH value accordingly. 
4)	Formulations may be saved for future use: By selecting the Save Formulation Database File menu option the formulation, including any added acids may be saved for future use. Saving the file will generate a new â€˜.matâ€™ database file which can be opened by the user. The formulation will be considered as one ingredient if reopened in IngredientDB. 
5)	Additional notes: 
  -	The buffer capacity curve is presented to help the user visualize the total buffering of a food ingredient mixture. The dark shaded area represents the tBeta value (see reference #1, above) for the mixture, the vertical lines represent individual buffers present in the mixture that contribute to the pH of the formulation. The dotted line and lightly shaded region represents the buffering of water, which by definition has a tBeta value of zero. 
  -	The Trim button on the graph will plot the ingredient formulation without buffering due to weakly acidic hydroxyl groups that are typically present on sugar molecules. This will also update total buffering and buffering percentages (but will not change formulation pH). For more information on the effects of sugars on pH and buffering, see reference #3, above. Note that the trim button does not adjust the buffering of added acid solutions, which were modeled from known pK values. 
  -	All ingredients used for a formulation must be titrated using the same NaCl% (based on ionic strength). An error message will appear if the user attempts to combine ingredients titrated with different salt concentrations.  
  -	The Clear button over the ingredient or acid table will remove all ingredients or change all acid concentrations to zero (respectively). 
  -	If there are no ingredients in the ingredient table the pH of acid solutions can be determined by changing the concentrations in the acid table, and the NaCl% box may be edited, changing the acid pK values accordingly (see reference 3 above for more details). However, once ingredients are added to the acid table the NaCl% will revert to the salt concentration used for the ingredient titrations. 
Contact information: 

For further help or questions please contact Dr. Fred Breidt: fred.breidt@usda.gov.   
