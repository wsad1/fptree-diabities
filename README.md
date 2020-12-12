# fptree-diabities
fp tree is a data mining algorithm to mine frequent item sets form data. This project implements fp tree on a diabetes database.

The dataset contains 9 features , some of the features are continuous or have to many values we need to discretize them before applying 
the fp tree algorithm, an entropy based clustering is used for this.

after getting frequent item sets we apply apriori algorithm to obtain rules.

References: 
Book
Ch 6 of Introduction to data Mining Pang,Seinbach,Kumar.
Ch 3 of Machine Learning, Tom M Mitchell
