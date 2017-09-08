#include <iostream>
#include <cstdlib>
#include <bits/stdc++.h>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#define pii pair<int,int>
using namespace std;
float minconf=0.7;
ofstream ft,fo;
struct node
{
	map<int,struct node*> m;
	int counter;
	int value;
	int flag;
	int pole;
	struct node *parent,*l;
};
struct tree
{
	node* nlinks[60];
	node *troot;
};

struct freq_item_set
{
    vector<int> v;
    int cnt;
};

vector<struct freq_item_set> arr[10];

int thresh_count =200;

int mdata[]={4,4,4,4,8,4,8,4};
string str[]={"pregnant1","pregnant2","pregnant3","pregnant4","PGC1","PGC2","PGC3","PGC4","BP1","BP2","BP3","BP4",
"Triceps1","Triceps2","Triceps3","insulin1","insulin2","insulin3","insulin4","insulin5","insulin6","insulin7","insulin8",
"BM1","BM2","BM3","BM4","Dia_pedigree_fun1","Dia_pedigree_fun2","Dia_pedigree_fun3","Dia_pedigree_fun4","Dia_pedigree_fun5","Dia_pedigree_fun6","Dia_pedigree_fun7","Dia_pedigree_fun8","Dia_pedigree_fun9","Dia_pedigree_fun10","Dia_pedigree_fun11","Dia_pedigree_fun12",
"Age1","Age2","Age3","Age4","no_diabities","diabities"};
string str2[46] ;
vector<int> data1[768];
float data[768][9];
int it_count[50];
void encode()
{
	//int data[495][17];
	std::ifstream file("diabities.txt");
	string line;
	int j,i;
	//cout<<"hi \n";
	for(j=0;j<768;j++)
	{
	    int count=0;
	 	std::getline(file,line);
		std::stringstream linestream(line);
		string item;
		for(i=0;i<9;i++)

		{

			std::getline(linestream,item,',');
			char *a = const_cast<char*>(item.c_str());
            if(i==0 || i==8)
                data[j][i]=float(atoi(a));
            else
                data[j][i]=float(atof(a));
            if(i!=0 && i!=8 && (data[j][i]==0.0 || data[j][i]==0))
            {
                data[j][i]=-500.0;
            }
        }
	}
}
vector<int> split(int lo,int hi,int num,int count,vector <pair<float,float> > arr)
{
    //printf("+++hi/n");
    int i,upcount1=0,upcount0=0,mid=-1,tot=hi-lo+1,fcount;
    float p=(float (count)/float (tot));
    float entro=-(p*log2f(p))-((1-p)*log2f(1-p)),tentro;
    vector<int> l,r;
    for(i=lo;i<hi;i++)
    {
        if(arr[i].second==1.0)
        {
            upcount1++;
        }
        p=(float(upcount1)/float(i-lo+1));
        tentro=-(p*log2f(p))-((1.0-p)*log2f(1-p));
        tentro=(float(i-lo+1)/float(tot))*tentro;
        p=(float(count-upcount1))/float(hi-i);
        tentro+=(float(hi-i)/float(tot))*(-(p*log2f(p))-((1.0-p)*log2f(1.0-p)));
        if(entro-tentro>0.0)
        {
            mid=i;
            entro=tentro;
            fcount=upcount1;
        }
    }
    if(num-2>0 && mid!=-1)
    {
        if(mid>lo+2 && fcount>0 && (mid-lo+1)-fcount>0)
            l=split(lo,mid,num-2,fcount,arr);
        if(hi>mid+2 && count-fcount>0 && hi-mid-(count-fcount)>0)
            r=split(mid+1,hi,num-2,count-fcount,arr);
    }
    if(mid==-1)
    {
        return l;
    }
    l.push_back(mid);
    l.insert(l.end(),r.begin(),r.end());
    return l;
}
void actual_split(vector<float> dpoints,int i,int num)
{
    int j,temp;
    int r=0;
    for(j=0;j<768;j++)
    {
        if(data[j][i]!=-500.0)
        {
            for(r=0;r<dpoints.size();r++)
            {
                if(data[j][i]<=dpoints[r])
                    break;
            }
            temp=r;
            temp++;
            data1[j].push_back(num-temp);
        }
    }
}

void revalue(int sz)
{
    int i=0,j;
    pair<int,int> parr[60];
    map<int,int> mp;
    for(j=0;j<768;j++)
    {
        for(i=0;i<data1[j].size();i++)
        {
            parr[-data1[j][i]].first++;
            parr[-data1[j][i]].second=data1[j][i];
        }

    }
    sort(parr+1,parr+sz+1);
    for(i=sz,j=1;i>=1;i--,j++)
    {
        it_count[j]=parr[i].first;
        mp[parr[i].second]=j;
        str2[j]=str[-parr[i].second-1];
    }
    for(j=0;j<768;j++)
    {
        for(i=0;i<data1[j].size();i++)
        {
            data1[j][i]=mp[data1[j][i]];
        }
    }
}
void preprocess()
{
    int i,j,r=0,count1=0,count2=0,sz=0;
    vector<int> spoints;
    vector<pair<float,float> > temp;
    vector <float> dpoints;
    encode();
    for(j=0;j<768;j++)
    {
        if(data[j][8]==1.0)
        {
            count1++;
        }
    }
    for(i=0;i<8;i++)
    {
        for(j=0;j<768;j++)
        {
            if(data[j][i]!=-500)
                temp.push_back(make_pair(data[j][i],data[j][8]));
        }
        sort(temp.begin(),temp.end());
        spoints=split(0,temp.size()-1,mdata[i],count1,temp);
        for(r=0;r<spoints.size();r++)
        {
                dpoints.push_back(temp[spoints[r]].first);
        }
        actual_split(dpoints,i,count2);
        count2=count2-spoints.size()-1;
        sz+=spoints.size()+1;
        spoints.clear();
        dpoints.clear();
        temp.clear();
    }
    dpoints.push_back(0);
    dpoints.push_back(1);
    actual_split(dpoints,8,count2);
    revalue(sz+2);
}

node* copyfnc(node* actroot,node* root,node* nlinks[])
{

	node *temp = root;
	node *nroot = new node();
	map<int,struct node*>::iterator it;
	nroot->parent = actroot;
	nroot->value = temp->value;
	nroot->flag=0;
	nroot->pole=0;
	nroot->counter = temp->counter;
	nroot->l=nlinks[nroot->value];
	nlinks[nroot->value]=nroot;
	for(it=temp->m.begin();it!=temp->m.end();it++)
	{
		nroot->m[it->first] = copyfnc(nroot,temp->m[it->first],nlinks);
	}
	return nroot;
}
int traverse(node* root,int item_count[],int k)
{
    vector<int> ra;
    node *temp = root;
    if(temp->value==k)
    {
        temp->m.clear();
        return temp->counter;
    }
    map<int,struct node*>::iterator it;
    temp->counter=0;
    for(it=temp->m.begin();it!=temp->m.end();it++)
    {

        if(it->second->flag==0)
        {
             ra.push_back(it->first);
        }
        else
        {
            temp->counter+= traverse(temp->m[it->first],item_count,k);
            item_count[temp->value] += temp->counter;
        }
    }
    for(int i=0;i<ra.size();i++)
    {
        temp->m.erase(ra[i]);
    }
    return temp->counter;
}
void traverse1(node *root,int item_count[])
{
    map<int,struct node*>::iterator it;
    for(it=root->m.begin();it!=root->m.end();it++)
    {
        int val = root->m[it->first]->value;
        if(item_count[val]<thresh_count)
            root->m[it->first]->pole = -1;
        traverse1(root->m[it->first],item_count);
    }
    return;
}
void make_sets(node *root,list<pii> l)
{
    int sz = l.size();
    freq_item_set titem;
    titem.cnt = l.front().first;
    list<pii>:: iterator it1;
    for(it1=l.begin();it1!=l.end();it1++)
    {
        titem.v.push_back(it1->second);
    }
    arr[sz].push_back(titem);
	tree t1;
	for(int i=0;i<60;i++)
		t1.nlinks[i]=NULL;
	t1.troot = copyfnc(NULL,root,t1.nlinks);
	int k = l.front().second;
	node *temp=t1.nlinks[k];
	while(temp!=NULL)
    {
        node *temp1 = temp;
        while(temp1->value!=0)
        {
            temp1->flag=1;
            temp1 = temp1->parent;
        }
        temp = temp->l;
    }
	int item_count[60];
    memset(item_count,0,sizeof(item_count));
    traverse(t1.troot,item_count,k);
    traverse1(t1.troot,item_count);

    for(int i=1;i<60;i++)
    {
        if(item_count[i]>thresh_count)
        {
            l.push_front({item_count[i],i});
            make_sets(t1.troot,l);
            l.pop_front();
        }
    }
}
void traverse3(node * rot,vector<int> tj)
{
    map<int,struct node*>::iterator it;
    int j=0;
    node *temp=rot;
    for(it=temp->m.begin();it!=temp->m.end();it++)
    {
        j++;
        tj.push_back(j);
        ft<<"(";
        for(int r=0;r<tj.size();r++)
        {
            ft<<tj[r]<<",";
        }
        ft<<"):";
        ft<<it->second->value<<"{"<<it->second->counter<<"}\n";
        traverse3(it->second,tj);
        tj.pop_back();
    }
}
void make_rule_i(vector <int> l1,vector<int> r1,int right,int i,int k)
{
    int j;
    vector<int> l2;
    float count;
    if(right>=k-1)
    {
        for(int r=0;r<arr[k-right].size();r++)
        {
            if(l1==arr[k-right][r].v)
            {
                count=float(arr[k-right][r].cnt);
                break;
            }

        }
        if(float(arr[k][i].cnt)/count>minconf)
        {
            fo<<"{";
            for(j=0;j<l1.size();j++)
            {
                fo<<str2[l1[j]]<<" ";
            }
            fo<<"}"<<count<<"=>{";
            for(j=0;j<r1.size();j++)
            {
                fo<<str2[r1[j]]<<" ";
            }
            fo<<"}"<<arr[k][i].cnt;
            fo<<"\n";
            return;
        }
        else
        {
                return;
        }
    }
    else
    {
        if(right>0)
        {
            for(int r=0;r<arr[k-right].size();r++)
            {
                if(l1==arr[k-right][r].v)
                {
                    count=arr[k-right][r].cnt;
                    break;
                }
            }
            if(float(arr[k][i].cnt)/count>minconf)
            {
                fo<<"{";
                for(j=0;j<l1.size();j++)
                {
                        fo<<str2[l1[j]]<<" ";
                }
                fo<<"}"<<count<<"=>{";
                for(j=0;j<r1.size();j++)
                {
                        fo<<str2[r1[j]]<<" ";
                }
                fo<<"}"<<arr[k][i].cnt;
                fo<<"\n";
            }
            else
            {
                return ;
            }
        }
        for(j=0;j<l1.size();j++)
        {
            for(int r=0;r<l1.size();r++)
            {
                if(r==j)
                    continue;
                else
                    l2.push_back(l1[r]);
            }
            if(right==0 || l1[j]>r1[r1.size()-1])
            {
                r1.push_back(l1[j]);
                make_rule_i(l2,r1,right+1,i,k);
                r1.pop_back();
            }
            l2.clear();
        }
    }
}
void make_all_rule(int k)
{
    int i=0,j=0;
    vector <int> l1;
    vector <int> r1;
    for(i=0;i<arr[k].size();i++)
    {
        l1=arr[k][i].v;
        make_rule_i(l1,r1,0,i,k);
    }
}
int main()
{
    fo.open("out2.txt");
    ft.open("out1.txt");
	int tot_items = 45;
	preprocess();
	tree t ;
	for(int i=0;i<60;i++)
		t.nlinks[i]=NULL;
	t.troot = new node();
	t.troot->value = 0;
	map<int,struct node*>::iterator it;
	for(int i=0;i<768;i++)
	{
		node *temp = t.troot;
		sort(data1[i].begin(),data1[i].end());
		for(int j=0;j<data1[i].size();j++)
		{
			if(it_count[data1[i][j]]>thresh_count)
			{
                int k = data1[i][j];
                it = temp->m.find(k);
                if(it==temp->m.end())
                {
                    temp->m[k] = new node();
                    temp->m[k]->parent = temp;
                    temp->m[k]->value = k;
                    temp->m[k]->counter = 1;
                    temp->m[k]->flag=0;
                    temp->m[k]->pole=0;
                    temp->l = t.nlinks[k];
                    t.nlinks[k] = temp;
                }
                else
                {
                    temp->m[k]->counter++;
                }
                temp=temp->m[k];
            }

        }
	}
	vector <int> tj;
	traverse3(t.troot,tj);
	for(int i = tot_items;i>0;i--)
	{
		if(it_count[i]>thresh_count)
		{
			list<pii> l1;
			l1.push_front({it_count[i],i});
			make_sets(t.troot,l1);
			l1.pop_front();
		}
	}
	for(int i=2;i<10;i++)
	{
	    if(arr[i].size()>0)
	    {
	        make_all_rule(i);
	    }
	}
	fo.close();
	ft.close();
    cin>>tot_items;
}
