# fptree-diabities
fp tree is a data mining algorithm to mine frequent item sets form data sets.the project implements fp tree on a database related to 
diabetes to get frequent item sets, we can experiment with different support counts.

The dataset contains 9 features , some of the features are continuous or have to many values we need to discretize them before applying 
the fp tree algorithm,i used entropy based clustering

after getting frequent item sets we apply apriori algorithm to obtain rules.

References: 
Book
Ch 6 of Introduction to data Mining Pang,Seinbach,Kumar.
Ch 3 of Machine Learning, Tom M Mitchell
