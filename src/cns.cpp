#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <string.h>
#include <algorithm>
#include <stdint.h>
#include <list>
#include <stack>
#include <limits.h>
#include <vector>
#define NINF INT_MIN
using namespace std;

// Default length, can override at runtime
int Kmer_Len = 30;
string CDG_Filename = "cdg.dot"; //for output of compressed de Bruijn graph in dot format

int DEBUG = 0;
int VERBOSE = 0;
int VERIFY = 0;

typedef uint32_t treeint;
typedef uint64_t treeintLarge;

treeintLarge skippedbases = 0;
treeintLarge skippedextensions = 0;

treeint numKmerLens = 0;
treeint numNodesWithTable = 0;
treeint numEntriesInAuxTables = 0;

// Settings for linear time alg
bool FORCEROOT = false;
bool DOJUMP = true;
bool DOINTERNALSKIP = true;
bool DOPHASETRICK = true;

bool MEM = false;

const int basecount = 7;
int b2i(char base){
    switch (base){
        case '$' : return 0;
        case 'A' : return 1;
        case 'C' : return 2;
        case 'G' : return 3;
        case 'N' : return 4;
        case 'T' : return 5;
        case '#' : return 6;

        default:
            cerr << "Unknown base: " << base << endl;
            return b2i('N');
    };
}

class MerVertex_t;

class SuffixNode{
    public:
        static treeintLarge s_nodecount;

    SuffixNode(
            treeint  s, //uint32_t
            treeint  e,
            SuffixNode * p,
            SuffixNode * x)
            : m_start(s),
              m_end(e),
              m_parent(p){//,
        s_nodecount++;
        m_isCopy = false;
        m_suffixTable = new SuffixNode*[1]; //for first suffix link, replaced with a larger table later

        m_suffixTable[0] = x;
        m_strdepth = 0;

        for (int i = 0; i < basecount; i++){
            m_children[i] = NULL;
        }
    }

    //copy constructor
    SuffixNode(const SuffixNode & node){
        m_start = node.m_start;
        m_end = node.m_end;
        m_SA_start = node.m_SA_start;
        m_SA_end = node.m_SA_end;

        m_parent = node.m_parent;
        for (int i = 0; i < basecount; i++){
            m_children[i] = NULL;
        }
        m_isCopy = true;

        m_suffixTable = node.m_suffixTable;
        m_LMAproximityTable= node.m_LMAproximityTable;
        m_LMAprox_nodeTable = node.m_LMAprox_nodeTable;

        m_MEM = node.m_MEM;
        m_LMA = node.m_LMA;
        m_strdepth = node.m_strdepth;
    }

    void createSuffixLinkTable(treeint numSuffixLinks, bool otherTables = true){
        if(numSuffixLinks==0)
            numSuffixLinks = 1; //min number of entries in table

        m_suffixTable = new SuffixNode*[numSuffixLinks];

        if(otherTables){
            numNodesWithTable++;

            m_LMAproximityTable = new treeint[numSuffixLinks];
            m_LMAprox_nodeTable = new SuffixNode*[numSuffixLinks];
        }else{
            m_LMAproximityTable = NULL;
            m_LMAprox_nodeTable = NULL;

        }
    }

    void deleteSuffixLinkTable(){
        SuffixNode * suffixLink = m_suffixTable[0];

        if(m_suffixTable){
            delete [] m_suffixTable;
        }
        if(m_LMAproximityTable){
            delete [] m_LMAproximityTable;
            m_LMAproximityTable = NULL;
        }
        if(m_LMAprox_nodeTable){
            delete [] m_LMAprox_nodeTable;
            m_LMAprox_nodeTable = NULL;

        }
        createSuffixLinkTable(1, false);
        m_suffixTable[0] = suffixLink;
    }

    ~SuffixNode(){
        if(!m_isCopy){
            if(m_suffixTable)
                delete [] m_suffixTable;
            if(m_LMAproximityTable)
                delete [] m_LMAproximityTable;
            if(m_LMAprox_nodeTable)
                delete [] m_LMAprox_nodeTable;
        }

        for (int i = 0; i < basecount; i++){
            if (m_children[i]) { delete m_children[i]; }
        }
    }

    string str(const string & s){
        return s.substr(m_start, m_end-m_start+1);
    }

    treeint len(int i=-1){
        if (i != -1){
            if (i < m_end){
                return i - m_start + 1;
            }
        }
        return m_end - m_start + 1;
    }

    ostream & printLabel(ostream & os, const string & str){
        string seq = str.substr(m_start, m_end-m_start+1);

        if (m_start == m_end && m_start == 0){
            os << "\"ROOT\"";
        }else{
            os << "\"" << seq;
            if(m_MEM){
                os << " MEM ";
            }
            os << ":"
               << " [" << m_start
               << "," << m_end << "]"
               << "}\"";
        }
        return os;
    }


    ostream & printNodeLabel(ostream & os){
        os << m_start << m_end << m_strdepth ;
        return os;
    }

    ostream & printEdgeLabel(ostream & os, const string & str){
        string seq = str.substr(m_start, m_end-m_start+1);
        os << "\"" << seq << (MEM?(m_MEM ? " MEM " :""):"")
           << " SA [" << m_SA_start
           << ", " << m_SA_end <<"]"
           << " [" << m_start
           << "," << m_end << "]\"";

        return os;
    }

    //returns lowest marked ancestor (if this node is marked, it is the lowest)
    //if no marked ancestor, NULL
    SuffixNode* LMA(){
        return m_LMA;
    }

    //return num suffix links in table based on strdepth
    inline int nodeNumSuffixLinks(){
        if(m_strdepth == 0)
            return 0; //for root

        return (int)(ceil(log(m_strdepth) / log(2)));
    }

    treeint  m_start;
    treeint  m_end;
    treeint  m_strdepth;

    bool m_isCopy; //set to true for copy constructor. then suffix link table will not be deleted when the function returns
    bool m_MEM;
    bool m_prevChar[basecount];

    typedef SuffixNode * SuffixNodePtr;

    SuffixNode * m_parent;
    SuffixNode * m_children [basecount];
    SuffixNodePtr * m_suffixTable;

    treeint * m_LMAproximityTable;  //[i] stores smallest m_LMAproximity encountered in m_suffixTable[i]
    SuffixNodePtr * m_LMAprox_nodeTable; //stores LMAnode corresponding to m_LMAproximityTable entries

    SuffixNode * m_LMA;
    //interval in suffix array that corresponds to this node
    treeint m_SA_start;
    treeint m_SA_end;
};
treeintLarge SuffixNode::s_nodecount(1);


ostream & operator<< (ostream & os, SuffixNode * n){
    return n->printNodeLabel(os);
}


class SuffixTree{
public:
    SuffixTree(const string & s)
            : m_nodecount(0), m_string(s), m_maxMEMstrdepth(0){
        m_root = new SuffixNode(0,0,NULL,NULL);
        m_root->m_suffixTable[0] = m_root;
        m_suffixArray = new treeint[s.length()];
    }

    treeint * m_suffixArray;
    SuffixNode * m_root;
    treeintLarge m_nodecount;
    string m_string;
    treeint m_maxMEMstrdepth; //max strdepth at any MEM node; only nodes whose strdepth \leq maxMEMstrdepth have suffix link tables and skipped MEM tables, each one is the size of their strdepth
    //std::forward_list<SuffixNodeMark> nodesWithSuffixSkips;
    SuffixNode ** nodesWithSuffixSkips;
    treeint nextNodeWithSuffixSkips; //next position in array to fill
    treeint numNodesWithSuffixSkips;


    void resetMaxMemStrdepth(){
        m_maxMEMstrdepth = 0;
    }

    void dumpNode(SuffixNode * node){
        int children = 0;
        for (int i = 0; i < basecount; i++){
            SuffixNode * child = node->m_children[i];
            if (child){
                children++;

                cout << " " << node << "->" << child;

                cout << " [minlen=" << child->len() << ", label=";
                child->printEdgeLabel(cout, m_string) << "]" << endl;

                dumpNode(child);
            }
        }

        if (node->m_suffixTable[0]){
            cout << " " << node << " -> " << node->m_suffixTable[0]
                 << " [style=dotted, constraint=false]" << endl;
        }

        if (children == 0){
            cout << " " << node << " [shape=box, label=";
            node->printLabel(cout, m_string) << "]" << endl;
        }else{
            cout << " " << node << " [label=";
            node->printLabel(cout, m_string) << "]" << endl;
        }
    }

    void dumpTree(){
        cerr << "Dumping tree" << endl;
        cout << "digraph G {" << endl;
        cout << " size=\"7.5,10\"" << endl;
        cout << " center=true" << endl;
        cout << " label=\"Suffix tree of \'" << m_string << "\' len:"
             << m_string.length()-1 << " nc:"
             << m_nodecount << "\"" << endl;
        dumpNode(m_root);
        cout << "}" << endl;
    }

    void dumpNodeText(ostream & out, SuffixNode * n, treeint depth){
        for (int b = 0; b < basecount; b++){
            if (n->m_children[b]){
                for (treeint i = 0; i < depth; i++){
                    out << " ";
                }
                out << " ";
                out << n->m_children[b]->str(m_string) << endl;
                dumpNodeText(out, n->m_children[b], depth+1);
            }
        }
    }

    void dumpTreeText(ostream & out){
        out << "Suffix Tree len=" << m_string.length()-1 << endl;
        out << "String: \"" << m_string << "\"" << endl;
        out << "+" << endl;
        dumpNodeText(out, m_root, 0);
    }

    void dumpTreeSorted(ostream & out, SuffixNode * node, const string & pathstring){
        int c = 0;

        string mystring = node->str(m_string);
        string ps(pathstring);
        ps.append(mystring);

        for (int i = 0; i < basecount; i++){
            if (node->m_children[i]){
                c++;
                dumpTreeSorted(out, node->m_children[i], ps);
            }
        }

        if (c == 0){
            out << ps << endl;
        }
    }


    SuffixNode * createNode(treeint s, treeint e, SuffixNode * p, SuffixNode * x){
        SuffixNode * retval = new SuffixNode(s, e, p, x);
        m_nodecount++;
        return retval;
    }

    //error to call this function with k<1.
    bool fillKhopSuffixLinkedList(SuffixNode * node, treeint k){
        bool hasSuffixSkip = false;
        //calculate num suffix links at node: log base 2 of strdepth at node
        int nodeSuffixLinks =  node->nodeNumSuffixLinks();

        int children = 0;
        for (int b = 0; b < basecount; b++)
            if(node->m_children[b])
                children++;

        //SM added last clause 8/6/14
        if(children == 0 || node->m_strdepth > m_maxMEMstrdepth || node->m_strdepth < Kmer_Len)  //no tables to fill for this node
            return hasSuffixSkip;

        if(nodeSuffixLinks > k){
            hasSuffixSkip = true;
            if(node->m_suffixTable[k-1]){
                //calculate suffix link that is 2^k hops away from node
                SuffixNode * kNeighbor = node->m_suffixTable[k-1]; //really the (k-1)Neighbor
                int kNeighborSuffixLinks =  kNeighbor->nodeNumSuffixLinks();

                if(kNeighborSuffixLinks > k-1){
                    node->m_suffixTable[k] = kNeighbor->m_suffixTable[k-1];

                    if(kNeighbor->m_LMAproximityTable[k-1] < node->m_LMAproximityTable[k-1]){
                        node->m_LMAproximityTable[k] = kNeighbor->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = kNeighbor->m_LMAprox_nodeTable[k-1];
                    }else{
                        node->m_LMAproximityTable[k] = node->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[k-1];
                    }
                }else{
                    node->m_suffixTable[k] = NULL;
                    node->m_LMAproximityTable[k] = node->m_LMAproximityTable[0];
                    node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[0];
                    //return; //chidren won't be any deeper as far as strdepth //SM added 6/12/14
                }
            }else{  //only applies to root
                node->m_suffixTable[k] = NULL;
            }
        }
        return hasSuffixSkip;
    }

    //error to call this function with k<1.
    void fillKhopSuffixNode(SuffixNode * node, treeint k, bool fillLinkedList){
        int nodeSuffixLinks =  node->nodeNumSuffixLinks();

        int children = 0;
        for (int b = 0; b < basecount; b++)
            if(node->m_children[b])
                children++;


        if(children == 0 || node->m_strdepth > m_maxMEMstrdepth )  //no tables to fill for leaf
            return;


        if(nodeSuffixLinks > k){
            if(fillLinkedList){
                //nodesWithSuffixSkips.push_front(node);  //SM added 6/16/14
                nodesWithSuffixSkips[nextNodeWithSuffixSkips++] = node;
            }

            numNodesWithSuffixSkips++;

            if(node->m_suffixTable[k-1]){
                //calculate suffix link that is 2^k hops away from node
                SuffixNode * kNeighbor = node->m_suffixTable[k-1]; //really the (k-1)Neighbor
                int kNeighborSuffixLinks =  kNeighbor->nodeNumSuffixLinks();

                if(kNeighborSuffixLinks > k-1){
                    node->m_suffixTable[k] = kNeighbor->m_suffixTable[k-1];

                    if(kNeighbor->m_LMAproximityTable[k-1] < node->m_LMAproximityTable[k-1]){
                        node->m_LMAproximityTable[k] = kNeighbor->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = kNeighbor->m_LMAprox_nodeTable[k-1];
                    }else{
                        node->m_LMAproximityTable[k] = node->m_LMAproximityTable[k-1];
                        node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[k-1];
                    }
                }else{
                    node->m_suffixTable[k] = NULL;
                    node->m_LMAproximityTable[k] = node->m_LMAproximityTable[0];
                    node->m_LMAprox_nodeTable[k] = node->m_LMAprox_nodeTable[0];
                }
            }else{  //only applies to root
                node->m_suffixTable[k] = NULL;
            }
        }
        for (int b = 0; b < basecount; b++){
            if(node->m_children[b]){
                fillKhopSuffixNode(node->m_children[b], k, fillLinkedList);
            }
        }
    }



    //recursively perform DFS to calc k-hop suffix links for each node
    void fillKhopSuffix(treeint k, bool fillLinkedList)
    {
        if(fillLinkedList) //dynamically allocate array
        {
            nodesWithSuffixSkips = new SuffixNode*[m_nodecount];
            nextNodeWithSuffixSkips = 0; //array position to fill next
        }

        fillKhopSuffixNode(m_root, k, fillLinkedList);
    }



    //this function uses linked list instead of multiple recursive DFSs over suffix tree
    //this function populates the table m_SuffixLinkTable and associated skipped LMA tables
    //entry i is the suffix link that is 2^i hops away
    //traverse suffix tree floor(log n) = m_numSuffixLinks times to populate the table
    void fillSuffixTableNonRecursively()
    {
        //timeval starttime;
        //timeval endtime;
        bool hasSuffixSkip;
        treeint nextSNodeSetup; //where up to in array
        treeint nextSNodeReplace; //where to copy to in array
        treeint totalSNodesWithSuffixSkips; //where working array ends
        //forward_list<SuffixNodeMark>::iterator it;

        treeint seconds;
        treeint microseconds;
        double elapsed;



        //begin with 1 since m_suffixLInkTable[0] is calcuated as part of Ukkonen's construction algorithm
        unsigned int numSuffixLinks = (int)( ceil(log(m_maxMEMstrdepth) / log(2))); //floor(log(m_maxMEMstrdepth - Kmer_Len) / log(2));

        cerr<<" numSuffixLinks="<<numSuffixLinks<<endl;


        //fill first suffix skip entries by traversing entire tree and filling nodesWithSuffixSkips linked list with only nodes that have suffix skips
        int k = 1;

        numNodesWithSuffixSkips = 0;

        //gettimeofday(&starttime, NULL);

        fillKhopSuffix(k, true);


        totalSNodesWithSuffixSkips = nextNodeWithSuffixSkips;

        cerr<<" total nodes with suffix skips ="<<nextNodeWithSuffixSkips<<endl;
        cerr<<" numNodesWithSuffixSkips = "<<numNodesWithSuffixSkips<<endl;
        cerr<<" now using linked lists"<<endl;


        for(k = 2; k < numSuffixLinks; k++)
        {
            numNodesWithSuffixSkips = 0;

            //gettimeofday(&starttime, NULL);

            nextSNodeReplace = 0;

            for( nextSNodeSetup = 0; nextSNodeSetup<totalSNodesWithSuffixSkips; nextSNodeSetup++)
            {
                hasSuffixSkip = fillKhopSuffixLinkedList(nodesWithSuffixSkips[nextSNodeSetup], k);
                if(hasSuffixSkip)
                {
                    nodesWithSuffixSkips[nextSNodeReplace++] = nodesWithSuffixSkips[nextSNodeSetup];
                    numNodesWithSuffixSkips++;
                }
                //else
                //nothing to do: loop advances nextSNodeSetup

            }

            totalSNodesWithSuffixSkips = numNodesWithSuffixSkips;



            cerr<<" numNodesWithSuffixSkips = "<<numNodesWithSuffixSkips<<endl;


        }

        delete [] nodesWithSuffixSkips;

    }


    void preprocessLMAnode(SuffixNode * node, SuffixNode * markedNode)
    {
        if(node->m_MEM && node->m_strdepth >= Kmer_Len)
        {
            markedNode = node;
        }
        node->m_LMA = markedNode;

        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                preprocessLMAnode(node->m_children[b], markedNode);
            }
        }


    }

    //mark internal nodes that represent MEMs of length >=minMEM by setting m_MEM = true
    void preprocessLMA( )
    {
        //recursively perform DFS on suffix tree and set LMA to point to lowest ancestor that is marked as a MEM
        preprocessLMAnode(m_root, NULL);
    }


    //recursively unmark all MEM nodes
    //can stop when reach a node that is long enough to be a MEM but isn't
    void unmarkMEMnode(SuffixNode * node, int minMEM)
    {
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
                unmarkMEMnode(node->m_children[b], minMEM);
        }

        node->m_MEM = false;
        node->m_LMA = NULL;

        node->deleteSuffixLinkTable();


    }

    void unmarkMEMnodes(int minMEM)
    {
        unmarkMEMnode(m_root, minMEM);  //root cannot be marked node
    }


    //recursively mark node if it's a MEM, fill suffix array along the way with same DFS
    void markMEMnode(SuffixNode * node, treeint minMEM, treeint * nextSAentry, bool setupSA)
    {
        int children = 0;
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                children++;
            }
        }

        if(setupSA)
        {
            //when visit a node the first time: set string depth in m_strdepth to reflect the length of the labels on the path from the root to this node's end
            if(node == m_root)  //set string depth = length of labels on path from root to this node's end
                node->m_strdepth = 0;
            else
                node->m_strdepth = node->m_parent->m_strdepth + node->len();

            //fill next suffix array element for this leaf
            if(children == 0)
            {
                m_suffixArray[*nextSAentry] = node->m_start - node->m_parent->m_strdepth;

                node->m_SA_start = *nextSAentry;
                node->m_SA_end = *nextSAentry;

                (*nextSAentry)++;
            }
        }

        //perform DFS recursively
        bool firstChild = true;


        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                markMEMnode(node->m_children[b], minMEM, nextSAentry, setupSA);
                for(int c = 0; c < basecount; c++)
                {
                    if(node->m_children[b]->m_prevChar[c])
                        node->m_prevChar[c] = true;

                }
                if(setupSA)
                {
                    if(firstChild) //first child
                    {
                        node->m_SA_start = node->m_children[b]->m_SA_start;
                        node->m_SA_end = node->m_children[b]->m_SA_end;
                        firstChild = false;
                    }
                    else
                    {
                        if(node->m_children[b]->m_SA_start < node->m_SA_start)
                            node->m_SA_start = node->m_children[b]->m_SA_start;

                        if(node->m_children[b]->m_SA_end > node->m_SA_end)
                            node->m_SA_end = node->m_children[b]->m_SA_end;
                    }
                }
            }
        }

        //when visit node the second time: set m_prevChar to true for its children's values (its value if leaf).  if more than one array element is set to true, set m_MEM to true if m_strdepth >= minMEM
        if(node == m_root)
        {
            node->m_MEM = false;
        }
        else if(children == 0) //leaf, cannot be MEM
        {
            node->m_MEM = false;
            treeint prevCharPos = node->m_end - node->m_strdepth;
            //if(prevCharPos >= 0) //ignore position 1, treat as $ before string so will be considered maximal (can't extend beyond beginning of string)
            {
                char prevChar;
                if(prevCharPos == 0)
                    prevChar = '$';
                else
                    prevChar = m_string[prevCharPos];

                int prevCharNum = b2i(prevChar);
                node->m_prevChar[prevCharNum] = true;
            }
        }
        else //internal node
        {
            if(node->m_strdepth >= minMEM) //see if this node is left maximal
            {

                int numPrevChars = 0;
                for(int b = 0; b < basecount; b++)
                {
                    if(node->m_prevChar[b])
                        numPrevChars++;
                }
                if(numPrevChars>=2) //can short-circuit OR if a child is a MEM
                {
                    node->m_MEM = true;

                    numKmerLens++;

                    if(node->m_strdepth > m_maxMEMstrdepth)
                        m_maxMEMstrdepth = node->m_strdepth;

                    //when visit MEM node the second time, also propagate up the start positions from the children, adjusted to account for the length of this node


                }
                else
                    node->m_MEM = false;
            }

        }

    }

    //mark internal nodes that represent MEMs of length >=minMEM by setting m_MEM = true
    //fill suffix array with same DFS in suffix tree
    void markMEMnodes(int minMEM, bool setupSA)
    {
        cerr<<" marking MEM nodes in suffix tree of at least "<<minMEM<<" bp"<<endl;
        treeint nextSAentry = 0;
        //recursively perform DFS on suffix tree and mark internal nodes that are left maximal as MEMs
        markMEMnode(m_root, minMEM, &nextSAentry, setupSA);
    }



    void createAuxTablesAtNode(SuffixNode * node)
    {
        SuffixNode * suffixLink = node->m_suffixTable[0];
        if(node->m_suffixTable != NULL)
            delete [] node->m_suffixTable;

        int children = 0;
        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
            {
                children++;
            }
        }

        if(node == m_root || children == 0 || node->m_strdepth > m_maxMEMstrdepth)  //leaves don't need tables for suffix links and skipped LMA
        {
            node->createSuffixLinkTable(1, false);
        }
        else
        {
            //don't go through all suffix links until root.  stop when strdepth of node is < Kmer_Len
            int numSuffixLinks = node->nodeNumSuffixLinks();
            if(numSuffixLinks==0)
                numSuffixLinks = 1; //min table size

            node->createSuffixLinkTable(numSuffixLinks, true);
        }

        node->m_suffixTable[0] = suffixLink;


        //from here comes from fillFirstLMAproxEntry_Node

        if(node != m_root && children > 0 && node->m_strdepth <= m_maxMEMstrdepth)
        {
            int minProx, thisProx, sProx;
            SuffixNode * snode = node->m_suffixTable[0];

            if(node->m_LMA)
                thisProx = node->m_strdepth - node->m_LMA->m_strdepth;  //this node's info is first entry
            else
                thisProx = node->m_strdepth;


            if(snode)
            {
                if(snode->m_LMA)
                {
                    sProx = snode->m_strdepth - snode->m_LMA->m_strdepth;  //this node's info is first entry
                }
                else
                {
                    sProx = snode->m_strdepth;
                }
                if(sProx < thisProx)
                    minProx = sProx;
                else
                    minProx = thisProx;
            }
            else
            {
                minProx = thisProx;
            }
            node->m_LMAproximityTable[0] = minProx;

            if(minProx == thisProx)
            {
                node->m_LMAprox_nodeTable[0] = node->LMA();
            }
            else
            {
                node->m_LMAprox_nodeTable[0] = snode->LMA();
            }
        } //to here

        for (int b = 0; b < basecount; b++)
        {
            if(node->m_children[b])
                createAuxTablesAtNode(node->m_children[b]);
        }
    }


    void createAuxTables()
    {
        numNodesWithTable = 0; //reset variable so can use in limited output as tables are created
        cerr<<" creating aux tables at nodes"<<endl;
        //recursively perform DFS on suffix tree and create suffix link and skipped LMA tables when appropriate
        createAuxTablesAtNode(m_root);
    }




    void buildUkkonen()
    {
        treeint len = m_string.length() - 1; // length of the string, not of the buffer
        char base = m_string[1];

        if (DEBUG)
        {
            cerr << "Building Ukkonen Tree for ";
            cerr << "string of len: " << len << endl;
        }

        // Construct T1
        SuffixNode * node = createNode(1, len, m_root, NULL);
        m_root->m_children[b2i(base)] = node;
        SuffixNode * firstleaf = node;
        SuffixNode * lastleaf = node;

        if (DEBUG)
        {
            cerr << "Phase 1 Child: ";
            node->printLabel(cerr, m_string) << endl;
        }

        treeint startj = 2;

        // phase i+1
        for (int i = 2; i <= len; i++)
        {
            DEBUG = 0;


            // Start at the last leaf created which will allow easy
            // access to the node for startj
            node = lastleaf;
            treeint nodewalk = 0;

            // Keep track of last internal nodes created in split so we can add suffix links
            SuffixNode * splitnode = NULL;

            if (!DOPHASETRICK)
            {
                startj = 2;
                node = firstleaf;
            }

            if (DEBUG)
            {
                char next = m_string[i];
                cerr << endl;
                cerr << i << ".0 " << "Phase " << i << " adding " << next << " starting with " << startj << endl;

                string beta = m_string.substr(1, i);
                cerr << i << ".1" << " Extension 1:  [implicit]" << endl;
            }

            for (treeint j = startj; j <= i; j++)
            {
                // Goal: Ensure S[j .. i] (beta) is in the suffix tree
                // Precondition: S[j-1 .. i] (alpha) is in the suffix tree "near" node
                //               All Internal nodes have a suffix link

                // Idea: 1) Remember where alpha is in the tree relative to node
                //       2) Walk up the tree w bases until we get to a node with a suffix link.
                //       3) Follow suffix link which shifts the path from S[j-1..i] to S[j..i]
                //       4) Walk down tree in new location ensuring S[i-w .. i] is in tree

                // Notes: 1) All internal nodes have a suffix link by next extension
                //        2) Any time we walk up to root, have to check S[j..i]
                //        3) Suffix [1..i] is always present so start extension j with 2

                treeint betapos = i; // The first position in string we need to check in tree

                if (DEBUG)
                {
                    cerr << endl;
                    string beta = m_string.substr(j, i-j+1);
                    cerr << i << "." << j << " Phase " << i << " Extension " << j << " bp:" << betapos << endl;

                    cerr << i << "." << j << "  Walking up from n:";
                    cerr << " nw: " << nodewalk << endl;
                }

                if (node == m_root)
                {
                    // If we are at root, we have to check the full string s[j..i] anyways
                }
                else
                {
                    if (nodewalk)
                    {
                        // partially walked down node->child, but didn't switch to child
                        // Match at i=6 on left... nodewalk=2, at 5 after suffix link
                        // 5 = i-2+1
                        //                 o ----- o
                        //               5 A       A 5  <-
                        //            -> 6 T       T 6

                        betapos -= nodewalk-1;

                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   Adjusted nw: " << nodewalk << endl;
                        }
                    }
                    else
                    {
                        // Exactly at a node or leaf.
                        // Walk up to parent, subtracting length of that edge
                        treeint len = node->len(i);
                        betapos -= len-1;
                        node = node->m_parent;

                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   Adjusted len: " << len << endl;
                        }
                    }

                    if (DEBUG)
                    {
                        cerr << i << "." << j << "   parent bp: " << betapos <<  " n:";
                        cerr<<endl;
                    }

                    if (node->m_suffixTable[0] == NULL)
                    {
                        // Subtract entire edge length
                        betapos -= node->len(i);
                        node = node->m_parent;

                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   grandparent bp: " << betapos << " n:";
                            cerr << endl;
                        }

                        if (node->m_suffixTable[0] == NULL)
                        {
                            cerr << "Missing suffix link!!! ";
                            exit(1);
                        }
                    }
                }

                // jump across suffix link
                node = node->m_suffixTable[0];
                if (node == m_root) { betapos = j; } // have to check full string

                if (DEBUG)
                {
                    cerr << i << "." << j << "  Starting to walk down from bp: " << betapos << " to " << i << " n:";
                    cerr << endl;
                }

                if (FORCEROOT && node != m_root)
                {
                    node = m_root;
                    betapos = j;

                    if (DEBUG)
                    {
                        cerr << i << "." << j << " AtRoot bp: " << betapos << endl;
                    }
                }

                bool done = false;
                startj = j+1; // assume this extension should be skipped in the next phase

                while ((betapos <= i) && !done)
                {
                    char base = m_string[betapos];
                    char b = b2i(base);
                    SuffixNode * child = node->m_children[b];

                    if (DEBUG)
                    {
                        cerr << i << "." << j << "  node betapos: " << betapos << "[" << base << "] n:";
                        cerr << endl;
                    }

                    if (!child)
                    {
                        if (splitnode && betapos == splitnode->m_start)
                        {
                            if (DEBUG)
                            {
                                cerr << i << "." << j << "   Add SL1: ";
                            }

                            splitnode->m_parent->m_suffixTable[0] = node;
                            splitnode = NULL;
                        }

                        SuffixNode * newnode = createNode(betapos, len, node, NULL);
                        node->m_children[b] = newnode;
                        lastleaf = newnode;

                        if (DEBUG)
                        {
                            cerr << i << "." << j << "   New Node: ";
                        }

                        node = newnode;

                        // This is the first base that differs, but the edgelength to
                        // i may be longer. Therefore set nodewalk to 0, so the entire
                        // edge is subtracted.
                        nodewalk = 0;
                        done = true;
                        break;
                    }
                    else
                    {
                        treeint nodepos = child->m_start;
                        nodewalk = 0;

                        char nodebase = m_string[nodepos];

                        if (nodebase != base)
                        {
                            char nb = m_string[nodepos];
                            cerr << "ERROR: first base on edge doesn't match edge label" << endl;
                            cerr << "       nb: " << nb << " base: " << base << endl;
                            exit(1);
                        }

                        // By construction, the string from j-1 to betapos to i-1
                        // must already by present in the suffix tree
                        // Therefore, we can skip checking every character, and zoom
                        // to exactly the right character, possibly skipping the entire edge

                        if (DOJUMP)
                        {
                            treeint mustmatch = i-1 - betapos + 1;
                            treeint childlen = child->len(i);

                            if (mustmatch >= childlen)
                            {
                                betapos += childlen;
                                nodepos += childlen;

                                skippedbases += childlen;

                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Edge Jump by: " << childlen << " new bp: " << betapos << " np: " << nodepos << endl;
                                }

                                if (nodepos != child->m_end+1)
                                {
                                    cerr << "ERROR: jump should have skipped entire edge, but didn't!" << endl;
                                    exit(1);
                                }
                            }
                            else if (mustmatch)
                            {
                                betapos += mustmatch;
                                nodepos += mustmatch;
                                nodewalk += mustmatch;

                                skippedbases += mustmatch;

                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Partial Jump by: " << mustmatch << " new bp: " << betapos << " np: " << nodepos << endl;
                                }

                                if (VERIFY)
                                {
                                    if (m_string[betapos-1] != m_string[nodepos-1])
                                    {
                                        cerr << "ERROR: jump should have matched at least the mustmatch-1 characters" << endl;
                                        cerr << "s[bp-1]: " << m_string[betapos-1] << " s[np-1]: " << m_string[nodepos-1] << endl;
                                        exit(1);
                                    }
                                }
                            }
                        }

                        while (nodepos <= child->m_end && betapos <= i)
                        {
                            nodebase = m_string[nodepos];

                            if (VERBOSE)
                            {
                                cerr << i << "." << j << "   child bp: " << betapos << "[" << m_string[betapos]
                                     << "] nb [" << nodebase << "]" << endl;
                            }

                            if (m_string[betapos] == nodebase)
                            {
                                if (splitnode && betapos == splitnode->m_start)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "   Add SL2: ";
                                    }

                                    splitnode->m_parent->m_suffixTable[0] = node;
                                    splitnode = NULL;
                                }

                                nodepos++; betapos++; nodewalk++;

                                if (betapos == i+1)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "    Internal edge match nw: " << nodewalk << endl;
                                    }

                                    if ((nodewalk == child->len(i)) && (child->m_end == len))
                                    {
                                        // we walked the whole edge to leaf, implicit rule I extension
                                        if (DEBUG)
                                        {
                                            cerr << i << "." << j << "    Leaf Node, Implicit Rule I Extension" << endl;
                                        }
                                    }
                                    else
                                    {
                                        // "Real" rule III implicit extension

                                        // The j-1 extension was the last explicit extension in this round
                                        // Start the next round at the last explicit extension
                                        if (DOPHASETRICK)
                                        {
                                            startj = j;

                                            treeint skip = startj - 2;

                                            if (DEBUG)
                                            {
                                                cerr << i << "." << j << "    Implicit Extension... start next phase at " << startj << ", saved " << skip << endl;
                                            }

                                            skippedextensions += skip;
                                        }

                                        if (DOINTERNALSKIP)
                                        {
                                            // Since we hit an internal match on a non-leaf, we know every other
                                            // extension in this phase will also hit an internal match.

                                            // Have to be careful since leafs get the full string immediately, but
                                            // they really have a Rule 1 extension

                                            treeint skip = i-j;

                                            if (DEBUG)
                                            {
                                                cerr << i << "." << j << "    Implicit Extension... skipping rest of phase, saved " << skip << endl;
                                            }

                                            skippedextensions += skip;
                                            j = i+1;
                                        }
                                    }

                                    done = true;
                                }
                            }
                            else
                            {
                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Spliting ";
                                }

                                // Split is a copy of the child with the end shifted
                                // Then adjust start of child
                                SuffixNode * split = createNode(child->m_start, nodepos-1, node, NULL);

                                split->m_children[b2i(nodebase)] = child;
                                child->m_start = nodepos;
                                child->m_parent = split;

                                if (DEBUG)
                                {
                                    cerr << " => ";
                                    //split->printLabel(cerr, m_string) << " + ";
                                    //child->printLabel(cerr, m_string) << endl;
                                }

                                node->m_children[b] = split;
                                node = split;

                                if (splitnode && betapos == splitnode->m_start)
                                {
                                    if (DEBUG)
                                    {
                                        cerr << i << "." << j << "   Add SL3: ";
                                        //splitnode->m_parent->printLabel(cerr, m_string) << " sl-> ";
                                        //node->printLabel(cerr, m_string) << endl;
                                    }

                                    splitnode->m_parent->m_suffixTable[0] = split;
                                    splitnode = NULL;
                                }

                                // Now create the new node
                                base = m_string[betapos];
                                b = b2i(base);
                                SuffixNode * newnode = createNode(betapos, len, split, NULL);
                                lastleaf = newnode;

                                split->m_children[b] = newnode;
                                splitnode = newnode;

                                node = newnode;

                                if (DEBUG)
                                {
                                    cerr << i << "." << j << "   Split New Node: ";
                                    //newnode->printLabel(cerr, m_string)
                                    cerr << endl;
                                }

                                // This is the first base that differs, but the edgelength to
                                // i may be longer. Therefore set nodewalk to 0, so the entire
                                // edge is subtracted.
                                nodewalk = 0;
                                done = true;
                                break;
                            }
                        }
                    }

                    if (!done) { node = child; }
                }
            }

            /*
             if (VERIFY)
             {
             verifySuffixLinks();
             }*/
        }
    }
};


SuffixTree * buildUkkonenSuffixTree(const string & s)
{
    SuffixTree * tree = new SuffixTree(s);
    tree->buildUkkonen();

    return tree;
}



void printHelp()
{

    cout << "findCNS [options] -file <fasta.file> -mem <CNS length>  -out <Output Prefix>"    << endl
         << "  -h               help"                                   << endl
         << "  -file <file>     Load sequence from file"                << endl
         << "  -mem <len>       Find CNSs at least this long "        << endl
         << "  -out         	   Prefix of the Output file " << endl;

    exit(0);
}



int OPT_DisplaySeq = 0;
int OPT_DisplayStarts = 1;
int OPT_SeqToDisplay = 8;
int OPT_DisplayStats = 0;
int OPT_PrintGraphs = 0;


// class for storing a mer and its neighbors
class MerVertex_t
{
    static treeintLarge NODECOUNT; // give all the nodes a unique id
public:
    MerVertex_t(treeint len)
            : node_m(NODECOUNT++) {
        length_m = len;
    }

    MerVertex_t(treeint startpos, treeint len)
            : node_m(NODECOUNT++) {
        starts_m.push_back(startpos);
        length_m = len;
    }

    MerVertex_t(treeint startpos, treeint len, bool noNodecount)
            : node_m(-1) {
        starts_m.push_back(startpos);
        length_m = len;
    }

    long   node_m; //unique node id
    vector<MerVertex_t *> successor_m;
    vector<MerVertex_t *> predecessor_m;
    vector<treeint> starts_m;
    treeint length_m;

    treeint addStartPos(treeint startpos)
    {
        starts_m.push_back(startpos);
        return starts_m.size()-1;
    }

    static treeintLarge getNodeCount()
    {
        return NODECOUNT;
    }

    treeint getNumEdges()
    {

        treeint numEdges = 0;
        numEdges = successor_m.size();
        return numEdges;
    }


    // Return the subsequence stored in the node
    string str(const string & seq)
    {
        string r = seq.substr(starts_m[0]-1, length_m);
        return r;
    }

    static void resetNodeCount()
    {
        NODECOUNT = 0;
    }
};

treeintLarge MerVertex_t::NODECOUNT=0;

//data type to wrap together a node pointer and which startpos index it represents
//this way the set of NodePos_t objects is searchable by startpos.  a vector of MerVertex_t is not searchable by startpos
class NodePos_t
{
public:
    MerVertex_t* nodePtr_m;
    treeint startPosIdx_m;

    NodePos_t(MerVertex_t* nptr, treeintLarge idx = 0){
        nodePtr_m = nptr;
        startPosIdx_m = idx;
    }

    void changeNodePtr(MerVertex_t* nptr)
    {
        nodePtr_m = nptr;
    }

    void changeStartIdx(int newIdx)
    {
        startPosIdx_m = newIdx;
    }

    bool operator ==(const NodePos_t * aNode)
    {
        if(getNodeBegin() == aNode->getNodeBegin())
            return true;
        else
            return false;
    }

    bool operator >(const NodePos_t * aNode)
    {
        if(getNodeBegin() > aNode->getNodeBegin())
            return true;
        else
            return false;
    }

    //get endPos based on startPosIdx of this nodePos
    treeint getNodeEnd() const
    {
        return getNodeBegin()+nodePtr_m->length_m-1;
    }

    //get startPos based on startPosIdx of this nodePos
    treeint getNodeBegin() const
    {
        return nodePtr_m->starts_m[startPosIdx_m];
    }

};


bool SortLens(const treeint & a, const treeint & b)
{
    return b < a;
}

bool SortStarts(const treeint & a, const treeint & b)
{
    return a < b;
}

bool SortNodePos(NodePos_t* a, NodePos_t* b)
{
    return a->nodePtr_m->starts_m[a->startPosIdx_m] < b->nodePtr_m->starts_m[b->startPosIdx_m];
}

bool SortMerVertex(MerVertex_t* a, MerVertex_t* b)
{
    return a->node_m < b->node_m;
}

typedef vector<NodePos_t *> NodeTable_t;
typedef vector<MerVertex_t *> RepeatNodeTable_t;
typedef vector<MerVertex_t *> UniqueNodeTable_t;



bool string_has_all_of_the_same_chars(const std::string& s) {
    //cout <<"Inside string_has_all_of...\n";
    return s.find_first_not_of(s[0]) == std::string::npos;
}
bool containsN(string str){
    return str.find("N") != std::string::npos;
}
bool containsX(string str){
    return str.find("X") != std::string::npos;
}
bool containsn(string str){
    return str.find("n") != std::string::npos;
}
bool containsx(string str){
    return str.find("x") != std::string::npos;
}

void seqPos(int *s, int y, int *n, int *z)
{
    int x=0, i=0;
    while(true)
    {
        *z = y - x;
        if((y-x) <= s[i]) break;
        // x = x + s[i] + 1;
        x = x + s[i];
        i++;
    }
    *n = i;
}


/* New printMEMnode() function */
void printMEMnode(SuffixNode * node, SuffixTree * ST, int *sLength, string seq[], string MEM[], int *index, int **matrix, int no_seq)
{
    //speedup: if !MEM and long enough for Kmer_Len, none of descendants are MEMs so don't bother
    if(!node->m_MEM && node->m_strdepth >= Kmer_Len )
        return;

    for (int b = 0; b < basecount; b++)
    {
        if(node->m_children[b])
        {
            printMEMnode(node->m_children[b], ST, sLength, seq, MEM, index, matrix, no_seq);
        }
    }

    int extendLeft = 0;

    if(node->m_MEM && node->m_strdepth >= Kmer_Len)
    {
        //todo:
        //this is a MEM, need to print it out
        if(node->m_parent != ST->m_root)
        {
            extendLeft = node->m_strdepth - node->len();
        }


        SuffixNode memNodeExt(*node);
        memNodeExt.m_start -= extendLeft;

        SuffixNode * memNode = &memNodeExt;
        memNode->m_strdepth = memNode->len();

        //first start
        //ST->m_suffixArray[memNode->m_SA_start], length
        //for each start in MEM node
        // cout<<"Detecting the MEMS......\n";
        string s = "";
        int t = *index;
        int d = t;

        int *count = (int *) malloc(sizeof(int) * no_seq);
        for(int i=0; i<no_seq; i++)
            count[i] = 0;

        for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
        {
            int y = ST->m_suffixArray[memNode->m_SA_start+i];
            int n=0, z=0;

            seqPos(sLength, y, &n, &z);
            //cout<<"y="<<y<<" n="<<n<<" z="<<z<<" depth"<<memNode->m_strdepth<<"\n";
            if (i==0)
            {

                s = seq[n].substr(z-1,memNode->m_strdepth);

            }
            if(!(string_has_all_of_the_same_chars(s) || containsn(s) || containsx(s) || containsN(s) || containsX(s))){
                count[n] +=1;
            }

        }

        int **ss = (int **) malloc(sizeof(int *) * no_seq);
        for(int i=0; i<no_seq; i++){
            ss[i] = (int *) malloc(sizeof(int) * count[i]+1);
        }

        int  *cc = (int *) malloc(sizeof(int) * no_seq);
        for(int i=0; i<no_seq; i++)
            cc[i] = 0;

        /*int  *X = (int *) malloc(sizeof(int) * no_seq);
        for(int i=0; i<no_seq; i++)
          X[i] = 0;*/

        for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
        {
            int y = ST->m_suffixArray[memNode->m_SA_start+i];
            int n=0, z=0;

            seqPos(sLength, y, &n, &z);
            //cout<<"y="<<y<<" n="<<n<<" z="<<z<<" depth"<<memNode->m_strdepth<<"\n";
            if (i==0)
            {
                s = seq[n].substr(z-1,memNode->m_strdepth);
            }


            if(!(containsN(s) || containsX(s) || containsN(s) || containsX(s))){
                int  *X = (int *) malloc(sizeof(int) * no_seq);
                for(int i=0; i<no_seq; i++)
                    X[i] = cc[i];
                int k = cc[n];
                int x = cc[0];
                ss[n][k] = z;
                //cc[n] = k + 1;
                X[n] = k+1;
                for(int i=0; i<no_seq; i++)
                    cc[i] = X[i];

            }

        }


        int n1=1, n2=no_seq;
        for(int i=0; i<no_seq; i++)
            if(cc[i] != 0) n1 *= cc[i];

        for(int i=0; i<n1; i++){
            int x = 0;
            //d++;
            for(int j=0; j<n2; j++){
                if(cc[j] == 0){
                    //cout<<"d="<<d<<"\n";
                    matrix[d][j] = 0;
                }
                else {
                    x = i%cc[j];
                    //cout<<"d="<<d<<"\n";
                    matrix[d][j] = ss[j][x];
                    x = i/cc[j];
                }
            }
            d++;
        }



        int aa = d - t;

        //if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s) || containsN(s) || containsX(s))){
        if(!(containsn(s) || containsx(s) || containsN(s) || containsX(s))){
            //cout<<"\n";
            for(int p=0; p<aa; p++){
                MEM[t]= s;
                t++;
                *index=t;
            }
        }

    }

}


/* New countMEMnode() function */
void countMEMnode(SuffixNode * node, SuffixTree * ST, int *sLength, string seq[],int *cnt, int no_seq)
{
    //cout<<"inside countMEM()\n";
    //speedup: if !MEM and long enough for Kmer_Len, none of descendants are MEMs so don't bother
    if(!node->m_MEM && node->m_strdepth >= Kmer_Len )
        return;

    for (int b = 0; b < basecount; b++)
    {
        //cout<<"\tinside for loop() with b= "<<b<<"\n";
        if(node->m_children[b])
        {
            //*cnt++;
            //int temp = *cnt;
            //cout <<"cnt:"<<temp<<"\n";
            countMEMnode(node->m_children[b], ST, sLength, seq, cnt, no_seq);
        }
    }
    int extendLeft = 0;
    if(node->m_MEM && node->m_strdepth >= Kmer_Len)
    {
        //todo:
        //this is a MEM, need to print it out
        if(node->m_parent != ST->m_root)
        {
            extendLeft = node->m_strdepth - node->len();
        }
        SuffixNode memNodeExt(*node);
        memNodeExt.m_start -= extendLeft;

        SuffixNode * memNode = &memNodeExt;
        memNode->m_strdepth = memNode->len();

        //first start
        //ST->m_suffixArray[memNode->m_SA_start], length
        //for each start in MEM node
        //cout<<"Detecting the MEMS......\n";
        //cnt++;
        int t = *cnt;



        int *count = (int *) malloc (sizeof(int)* no_seq);
        for(int i=0; i< no_seq; i++) count[i]=0;

        string s;
        for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
        {
            int y = ST->m_suffixArray[memNode->m_SA_start+i];
            int n=0, z=0;
            int xyz = 0;

            seqPos(sLength, y, &n, &z);
            //cout << "y: "<<y << " n: "<<n <<" z: "<< z<<"length:"<<memNode->m_strdepth<<"\n";
            if (i==0)
            {

                s = seq[n].substr(z-1,memNode->m_strdepth);
                //cout<<"s="<<s<<"\n";

            }
            //cout << "y: "<<y << " n: "<<n <<" z: "<< z<<"\n";
            //cout<< seq[n].substr(z-1,memNode->m_strdepth)<< " Seq_"<<n+1<<"\t"<<z<<"\t"<< memNode->m_strdepth << endl;
            //if(!(string_has_all_of_the_same_chars(s) || containsN(s) || containsX(s) || containsN(s) || containsX(s))){
            if(!(containsn(s) || containsx(s) || containsN(s) || containsX(s))){
                count[n]+=1;
            }



        }
        int flag = 0;
        for(int i=0; i<no_seq; i++)
            if(count[i] !=0 ) flag = 1;

        int k = 1;
        if(flag != 0) {
            //int k = 1;
            for(int i=0; i<no_seq; i++){
                if(count[i] == 0) count[i] += 1;
                k *= count[i];
            }
        }
        t += k;
        *cnt = t;



    }

}

void printMEMs(SuffixTree * tree, int *sLength, string seq[], string MEM[], int **matrix, int no_seq)
{
    //recursively perform DFS on suffix tree to find all MEMs
    int index=0;
    printMEMnode(tree->m_root, tree, sLength, seq, MEM, &index, matrix, no_seq);
    /*int cnt=0;
    countMEMnode(tree->m_root, tree, &cnt);
    cout<<"The no. of MEM nodes is:"<<cnt<<"\n";*/

}
// find single source longest distances in a DAG


// Graph is represented using adjacency list. Every node of adjacency list
// contains vertex number of the vertex to which edge connects. It also
// contains weight of the edge
class AdjListNode
{
    int v;
    int weight;
public:
    AdjListNode(int _v, int _w)  { v = _v;  weight = _w;}

    int getV()       {  return v;  }
    int getWeight()  {  return weight; }
};

class revAdjListNode
{
    int v;
    int weight;
public:
    revAdjListNode(int _v, int _w)  { v = _v;  weight = _w;}

    int getV()       {  return v;  }
    int getWeight()  {  return weight; }
};

// Class to represent a graph using adjacency list representation
class Graph
{
    int V;    // No. of verticesí

    // Pointer to an array containing adjacency lists
    list<AdjListNode> *adj;
    list<revAdjListNode> *revAdj;

    // A function used by longestPath
    void topologicalSortUtil(int v, bool visited[], stack<int> &Stack);
public:
    Graph(int V);   // Constructor

    // function to add an edge to graph
    void addEdge(int u, int v, int weight);

    // Finds longest distances from given source vertex
    void longestPath(int s, int **matrix, string *new_mems, int count, int no_seq, vector<int *> &LP_Matrix, vector<string> &LP_MEM);
};

Graph::Graph(int V) // Constructor
{
    this->V = V;

    adj = new list<AdjListNode>[V];
    revAdj = new list<revAdjListNode>[V];
}

void Graph::addEdge(int u, int v, int weight)
{
    AdjListNode node(v, weight);
    //if(u==35) cout<<"v ="<<v<<endl;
    revAdjListNode node1(u, weight);
    adj[u].push_back(node); // Add v to uís list
    revAdj[v].push_back(node1);

}

// A recursive function used by longestPath.
void Graph::topologicalSortUtil(int v, bool visited[], stack<int> &Stack)
{
    // Mark the current node as visited
    visited[v] = true;

    // Recur for all the vertices adjacent to this vertex
    list<AdjListNode>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {
        AdjListNode node = *i;
        if (!visited[node.getV()])
            topologicalSortUtil(node.getV(), visited, Stack);
    }

    // Push current vertex to stack which stores topological sort
    Stack.push(v);
}

// The function to find longest distances from a given vertex. It uses
// recursive topologicalSortUtil() to get topological sorting.
void Graph::longestPath(int s, int **matrix, string *new_mems, int count, int no_seq, vector<int *> &LP_Matrix, vector<string> &LP_MEM)
{

    int u = count;
    //cout<< "u = "<<u<<endl;
    list<AdjListNode>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i){
        //cout<< "I am here..."<<endl;
        //cout<<"i->getV()="<<i->getV()<<" i->getWeight()="<<i->getWeight()<<endl;

    }
    //exit(0);

    cout<<"Inside longest path .."<<endl;

    int ** mat = (int **) malloc(sizeof(int *) * count);
    for(int k=0; k <count; k++)
        mat[k] = (int *) malloc(sizeof(int)*no_seq + 1);
    for(int i=0; i<count; i++){
        int diff = 0;
        for(int j=0; j<no_seq; j++){
            mat[i][j] = matrix[i][j];
            for(int k=j+1; k<no_seq; k++){
                diff += fabs(matrix[i][j] - matrix[i][k]);
            }
            //mat[k][ind] = fabs(mat[k][2] - mat[k][3]) + fabs(mat[k][2] - mat[k][4]) + fabs(mat[k][3] - mat[k][4]);
        }
        mat[i][no_seq] = diff;
        //k++;
    }

    for(int k=0; k<count; k++){
        for(int j=0; j<no_seq+1; j++) cout << mat[k][j]<<" ";
        cout<<endl;

    }
    //exit(0);
    stack<int> Stack;
    int dist[V];

    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;

    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);
    //cout<<"for testing -- 2\n";
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++)
        dist[i] = NINF;
    dist[s] = 0;

    // Process vertices in topological order
    //cout<< "testing..while loop"<<endl;
    while (Stack.empty() == false)
    {
        // Get the next vertex from topological order
        int u = Stack.top();
        //cout<<"u="<<u<<endl;
        //cout<<"u->getWeight() "<<u->getWeight()<<endl;
        Stack.pop();

        // Update distances of all adjacent vertices
        list<AdjListNode>::iterator i;
        if (dist[u] != NINF)
        {
            //cout<<"dist[u] != NINF"<<u<<endl;
            for (i = adj[u].begin(); i != adj[u].end(); ++i){
                //cout<< "I am here..."<<endl;
                //cout<<"i->getV()="<<i->getV()<< "dist[u]="<<dist[u]<<" i->getWeight()="<<i->getWeight()<<endl;
                if (dist[i->getV()] < dist[u] + i->getWeight())
                    dist[i->getV()] = dist[u] + i->getWeight();

            }
        }
    }
    //cout<<"testing ends...while loop"<<endl;
    // Print the calculated longest distances
    int max = s;
    /*vector<int *> LP_Matrix;
    vector<string> LP_MEM;*/
    //FILE *MEMList = NULL;
    //char temp[100];
    for (int i = 0; i < V; i++)
    {
        //cout<<i<<"--"<<dist[i]<<endl;
        //(dist[i] == NINF)? cout << "INF ": cout << dist[i] << " ";
        if(dist[i] >= dist[s]) max = i;
    }
    //cout<<"The destination vertex is "<<max<<" with score "<<dist[max]<<"\n";

    int flag = 0;
    int start = max;
    int st = max;
    //cout<<max<<"("<<dist[max]<<")";
    list<revAdjListNode>::iterator i1;
    int last = 0;

    while(start != s){
        //cout << "Start:"<<start<<endl;
        int d = -999;
        for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
            if(dist[i1->getV()] > d){
                d = dist[i1->getV()];
                st = i1->getV() ;
            }
        }
        if(st == s) break;
        //cout<<"TESTING HERE "<<endl;
        //cout<<" <-( ";
        int z=0; int dd = 0; int st1 = 0;
        for(i1 = revAdj[start].begin(); i1 != revAdj[start].end(); ++i1){
            // cout<< " t_s: "<<i1->getV()<<"-"<<dist[i1->getV()]<<" t_e ";
            if(dist[i1->getV()] == d) {
                if(last < 1){
                    if(z == 0){
                        dd = mat[i1->getV()][no_seq];
                        st1 = i1->getV();
                    }
                    else {
                        if(mat[i1->getV()][no_seq] < dd){
                            dd = mat[i1->getV()][no_seq];
                            st1 = i1->getV();
                        }
                    }
                }
                //cout<< "TESTING HERE2 " << endl;
                /* added 10_14 */
                if(last >= 1){
                    if(z == 0){
                        //dd = mat[i1->getV()+1][5];
                        for(int i=0 ; i<no_seq; i++){
                            for(int j=i+1; j<no_seq; j++){
                                //cout<<"start = "<<start<<" i = "<<i<<" j = "<<j<<" i1->getV() = "<<i1->getV()<<endl;
                                int d1 = mat[start][i] - mat[i1->getV()][i];
                                int d2 = mat[start][j] - mat[i1->getV()][j];
                                dd += fabs(d1 - d2);
                                //cout<<" dd = "<<dd<<endl;
                            }
                        }
                        st1 = i1->getV();
                        //cout<<"TESTING HERE3 "<<endl;
                    }
                        //cout<<"TESTING HERE3 "<<endl;
                    else {
                        int ddd;
                        for(int i=0 ; i<no_seq; i++){
                            for(int j=i+1; j<no_seq; j++){
                                int d1 = mat[start][i] - mat[i1->getV()][i];
                                int d2 = mat[start][j] - mat[i1->getV()][j];
                                ddd += fabs(d1 - d2);
                            }
                        }
                        if(ddd< dd){
                            dd = ddd;
                            st1 = i1->getV();
                        }
                    }
                }
                z++;
            }
        }
        //cout<<"--st1:"<<st1<<endl;

        LP_Matrix.push_back(mat[st1]);
        LP_MEM.push_back(new_mems[st1]);
        start = st1;
        last++;
    }
}


void LP(int **matrix, string *new_mems, int count, int no_seq, vector<int *> &LP_Matrix, vector<string> &LP_MEM, string *seqName)
{

    cout << "Inside LP "<<count<<endl;
    int *ml = new int[count];
    for(int i = 0; i < count; i++){
        ml[i] = new_mems[i].length();
    }

    Graph g(count + 2);
    int newCount = count + 2;
    int ** adjMat = (int **) malloc(sizeof(int*)* newCount + 1);
    for(int i = 0; i < newCount; i++){
        adjMat[i] = (int *) malloc(sizeof(int) * 2 + 1);
        for(int j = 0; j < 2; j++) adjMat[i][j] = 0;
    }

    /*int **deg = (int **) malloc(sizeof(int *) * newCount + 1);
    for(int i = 0; i < newCount; i++){
      deg[i] = (int *) malloc(sizeof(int)*2 + 1);
      for(int j = 0; j < 2; j++) deg[i][j] = 0;
    } */




    int v1 = count, v2 = count + 1;
    int c3 = 0, s = 0;
    cout << "count" << count << endl;
    for(int p = 0; p < count; p++){

        //v1 -------> p Edge
        if(adjMat[p][0] == 0){
            g.addEdge(v1, p, ml[p]);
            adjMat[p][0] = 1;
        }

        // p -----> v2 Edge
        if(adjMat[p][1] == 0){
            g.addEdge(p, v2, 1);
            adjMat[p][1] = 1;
        }

        for(int q = p + 1; q < count; q++){
            int flag1 = 0, flag2 = 0, flag3 = 0;
            for(int r = 0; r < no_seq; r++){
                if(matrix[p][r] > matrix[q][r] && matrix[p][r] > (matrix[q][r] + ml[q])) flag1++;
                if(matrix[p][r] < matrix[q][r] && (matrix[p][r] + ml[p]) < matrix[q][r]) flag2++;
            }

            if(flag1 == no_seq ){
                // q ----> p Edge
                g.addEdge(q, p, ml[p]);

                // v1 -----> q Edge
                if(adjMat[q][0] == 0){
                    g.addEdge(v1, q, ml[q]);
                    adjMat[q][0] = 1;
                }
            }
            if(flag2 == no_seq){
                // Edge  p ----> q
                g.addEdge(p, q, ml[q]);

                // Edge q -----> v2
                if(adjMat[q][1] == 0){
                    g.addEdge(q, v2, 1);
                    adjMat[q][1] = 1;
                }
            }
        }
    }

    cout<< "Calling longest path ..."<<endl;
    g.longestPath(count, matrix, new_mems, count, no_seq, LP_Matrix, LP_MEM);
    cout<<"after calling longest path..."<<endl;
    delete[] adjMat;
    //delete[] deg;
    delete[] ml;

}



class Node
{
public:
    void AddLink(int id)
    {
        next.push_back(id);
    }

public:
    vector <int> next;
};

void FindAllPathsAt(vector <Node> &all_nodes, int id, vector < vector<int> > &all_paths, vector <int> tmp)
{
    tmp.push_back(id);

    if(all_nodes[id].next.size() == 0) {
        all_paths.push_back(tmp);
        return;
    }

    for(size_t i=0; i < all_nodes[id].next.size(); i++) {
        vector <int> tmp2(tmp);
        FindAllPathsAt(all_nodes, all_nodes[id].next[i], all_paths, tmp2);
    }
}

void PrintPaths(const vector < vector<int> > &all_paths, int **matrix, int count)
{
    int i,j;
    /*for(i=0; i<count-1; i++){
      for(j=0; j<5; j++)
        printf("%d ", matrix[i][j]);
      printf("\n");  
    }*/

    int num_rows = sizeof(matrix) / sizeof(matrix[0]);
    int num_cols = sizeof(matrix[0]) / sizeof(matrix[0][0]);

    //cout <<" rows: "<< num_rows << "cols:"<< num_cols; 
    ofstream myfile ("allPaths.txt");
    size_t max = 0;
    for(size_t i=0; i < all_paths.size(); i++) {
        if(max < all_paths[i].size()) max = all_paths[i].size();
    }
    //cout << "max=" << max;
    if (myfile.is_open()){

        for(size_t i=0; i < all_paths.size(); i++) {
            // Don't print node if it points to nothing
            if(all_paths[i].size() == 1) {
                continue;
            }

            if(all_paths[i].size() == max)
            {
                myfile<<all_paths[i][0];
                //cout << all_paths[i][0];

                int score1 = 0;
                int score2 = 0;
                int prev = 0;
                for(size_t j=1; j < all_paths[i].size(); j++) {
                    //if(all_paths[i][j].size() == max)
                    //cout << " -- > " << all_paths[i][j];
                    myfile<<"-"<<all_paths[i][j];
                    //if(j != all_paths[i].size()) score1 += matrix[all_paths[i][j]][0];
                    /*if(j != 1 || j != all_paths[i].size()){
                      int x = matrix[all_paths[i][j]][1] - matrix[prev][1];
                      int y = matrix[all_paths[i][j]][2] - matrix[prev][2];
                      int z = matrix[all_paths[i][j]][3] - matrix[prev][3];
                      score2 += fabs(x-y)+fabs(x-z)+fabs(y-z);

                    }*/
                    prev = all_paths[i][j];
                }
                myfile<< " ("<<score1<<")";
                myfile<< " ("<<score2<<")";
                myfile<<"\n";
                //cout << endl;
            }

            //cout << endl;
        }
    }
    myfile.close();
}

void allPath(char *fname, char *cID)
{

    FILE *instream = fopen(fname, "r");
    if(instream == NULL) {
        fprintf(stderr, "Unable to open file: %s\n", fname);
        exit(1);
    }

    printf("file name is %s\n", fname);
    char temp[100];
    int no_of_lines=0;
    while(fgets(temp, 100, instream)!= NULL) {
        no_of_lines++;
    }
    fclose(instream);
    int count=0;
    int **matrix = (int **) malloc(sizeof(int *)*no_of_lines +1);
    int i = 0;
    for(i=0; i<no_of_lines; i++)
        matrix[i] = (int *) malloc(sizeof(int)*5);
    //int matrix[100][4];
    int *ml = (int *) malloc(sizeof(int)*no_of_lines+1);
    //int ml[100];
    instream = fopen(fname, "r");
    while(fgets(temp, 100, instream)!= NULL) {

        temp[strlen(temp)-1]='\0';
        //printf("%s\n", temp);
        if(count!=0){
            int m = atoi(strtok(temp, "\t"));
            //	printf("%d\n", m);
            int l = atoi(strtok(NULL, "\t"));
            //printf("%d\n", l);
            int x = atoi(strtok(NULL, "\t"));
            //printf("%d\n", x);
            int y = atoi(strtok(NULL, "\t"));
            //printf("%d\n", y);
            //if(count==1) printf("%d %d %d %d", m,l,x,y);
            int z = atoi(strtok(NULL, "\t"));
            //printf("%d\n", z);
            //if(count==1)printf("%d %d %d %d %d", m,l,x,y,z);
            ml[count-1]=l;
            matrix[count-1][0]=l;
            matrix[count-1][1]=x;
            matrix[count-1][2]=y;
            matrix[count-1][3]=z;
        }
        count++;
    }
    printf("count=%d\n", count);
    FILE *f;
    char fn[] = "adj.txt";
    printf("fn =%s\n", fn);
    f=fopen(fn, "w");
    count--;
    //Graph g(count+2);
    vector <Node> all_nodes(count+2);
    int j=0;
    i = 0;
    int newCount = count+2;
    //int **adjMat = (int **) malloc(sizeof(int *)*count+1);
    int **adjMat = (int **) malloc(sizeof(int *)*newCount+1);
    for(i=0; i<newCount; i++){
        adjMat[i] = (int *) malloc(sizeof(int)*newCount+1);
        for(j=0; j<newCount; j++) adjMat[i][j]=0;
    }
    int **deg = (int **) malloc(sizeof(int *)*newCount+1);
    for(i=0; i<newCount; i++)
        deg[i] = (int *) malloc(sizeof(int)*2+1);
    for(i=0; i<newCount; i++){
        deg[i][0]=0;
        deg[i][1]=0;
    }
    int p,q,r;
    //cout << "cnt = "<<count<<"\n";
    int v1 = count, v2 = count+1;
    int c3=0, s=0;
    for(p=0; p<count; p++){
        //if(p<2) {printf("TESTING MATRIX %d %d %d\n", matrix[p][0], matrix[p][1], matrix[p][2]);}
        for(q=p+1; q<count; q++){
            int c2=0;
            int flag1 = 0, flag2=0, flag3=0;
            for(r=0; r<3; r++){
                if(matrix[p][r] > matrix[q][r]) flag1++;
                if(matrix[p][r] < matrix[q][r]) flag2++;
            }
            c3++;
            if(flag1 == 3 ){
                //if(p==0) printf("TESTING p = 0, q = %d\n", q);
                // q ----> p Edge
                //g.addEdge(q, p, ml[p]);
                all_nodes[q].AddLink(p);
                adjMat[q][p]=1;
                deg[q][1]+=1;  //out-degree of q
                deg[p][0]+=1;  // in-degree of p

                // v1 -----> q Edge
                if(adjMat[v1][q]==0){
                    //g.addEdge(v1, q, ml[q]);
                    all_nodes[v1].AddLink(q);
                    deg[q][0]+=1;
                    deg[v1][1]+=1;
                    adjMat[v1][q]=1;
                }

                // p -----> v2 Edge
                if(adjMat[p][v2]==0){
                    //g.addEdge(p, v2, 1);
                    all_nodes[p].AddLink(v2);
                    deg[v2][0]+=1;
                    deg[p][1]+=1;
                    adjMat[p][v2]=1;
                }
            }
            if(flag2 == 3){
                //if(p==0 ) printf("TESTING p = 0, q = %d\n", q);
                // Edge  p ----> q
                //g.addEdge(p, q, ml[q]);
                all_nodes[p].AddLink(q);
                adjMat[p][q]=1;
                deg[p][1]+=1;
                deg[q][0]+=1;

                // Edge v1 -----> p
                if(adjMat[v1][p] == 0){
                    //g.addEdge(v1, p, ml[p]);
                    all_nodes[v1].AddLink(p);
                    deg[p][0]+=1;
                    deg[v1][1]+=1;
                    adjMat[v1][p]=1;
                }


                // Edge q -----> v2
                if(adjMat[q][v2]==0){
                    //g.addEdge(q, v2, 1);
                    all_nodes[q].AddLink(v2);
                    deg[q][1]+=1;
                    deg[v2][0]+=1;
                    adjMat[q][v2]=1;
                }

            }
            //if(c3==1) s = p;
        }

    }

    /*for(i=0; i<newCount; i++){
      for(j=0; j<newCount; j++)
        fprintf(f, "%d ", adjMat[i][j]);
      fprintf(f, "\n");
    }  */


    /* printf("Printing in-degree and out-degrees\n");
     for(i=0; i<newCount; i++){
       printf("%d\t%d\t%d\n", i, deg[i][0], deg[i][1]);
       //printf("\n");
     }*/
    //g.longestPath(count, matrix, ml, cID);
    vector <int> tmp;
    vector <vector<int> > all_paths;
    FindAllPathsAt(all_nodes, v1, all_paths, tmp);
    cout<< "All paths at node "<< v1 << endl;
    PrintPaths(all_paths, matrix, count);
    /*for(i=0; i<count; i++){
     
        cout << "Following are longest distances from source vertex " << i <<" \n";
         g.longestPath(i);
         printf("\n");


    }*/
    //return 0;
}

/* checks if the string contains '#' (returns 1)*/
bool containsD(string str){
    return str.find("#") != std::string::npos;
}


void cartesian(vector<int> *memPos, int no_seq, vector<int *> &memTable, vector<string> &mems, string mem)
{
    int n1 = 1, n2 = no_seq;
    //int *arr = new int[no_seq];
    for (int i = 0; i < no_seq; i++)
        if (!memPos[i].empty()) n1 *= memPos[i].size();

    for (int i = 0; i < n1; i++) {
        int *arr = new int[no_seq];
        int x = 0;
        int operand = i;
        for (int j = 0; j < n2; j++) {
            if (memPos[j].empty()) {
                arr[j] = 0;
            }
            else {
                x = operand % memPos[j].size();
                //cout<<"at:"<< j<<": "<<memPos[j].at(x)<<endl;
                arr[j] = memPos[j].at(x) + 1;
                operand /= memPos[j].size();
            }
        }
        /*for (int i = 0; i < no_seq; i++) cout<< arr[i] << "\t";
        cout<<endl;*/
        memTable.push_back(arr);
        mems.push_back(mem);

    }
}

void printMEMnode_original(SuffixNode * node, SuffixTree * ST, string S, int *sLength, int numSeq, vector<int *> &memTable, vector<string> &mems)
{
    //speedup: if !MEM and long enough for Kmer_Len, none of descendants are MEMs so don't bother
    if(!node->m_MEM && node->m_strdepth >= Kmer_Len )
        return;

    for (int b = 0; b < basecount; b++)
    {
        if(node->m_children[b])
        {
            printMEMnode_original(node->m_children[b], ST, S, sLength, numSeq, memTable, mems);
        }
    }

    int extendLeft = 0;
    /* Added on 02/15/2015. starts here.... */
    int flag = 0;
    int count = 0;
    int sum = 0;
    for(int b = 0; b < basecount; b++)
    {
        if(node->m_prevChar[b]){
            sum = sum + b*pow(10, count);
            count++;
            //numPrevChars++;
        }
    }
    //cout<<"children "<<sum<<endl;
    for (int p = 0; p < basecount; p++)
    {
        int count1 = 0;
        int sum1 = 0;
        for(int b = 0; b < basecount; b++)
        {
            if(node->m_children[p] != NULL){
                if(node->m_children[p]->m_prevChar[b]){
                    sum1 = sum1 + b*pow(10, count1);
                    count1++;
                }
            }
        }

        if(sum1>10 && sum > 10 && sum1 == sum) {
            flag = 1;
            //cout<<"children1= "<<sum1<<endl;
        }

    }
    /* ends here....*/

    //if(node->m_MEM && node->m_strdepth >= Kmer_Len)
    if(node->m_MEM && node->m_strdepth >= Kmer_Len && flag == 0)
    {
        //todo:
        //this is a MEM, need to print it out
        if(node->m_parent != ST->m_root)
        {
            extendLeft = node->m_strdepth - node->len();
        }


        SuffixNode memNodeExt(*node);
        memNodeExt.m_start -= extendLeft;

        SuffixNode * memNode = &memNodeExt;
        memNode->m_strdepth = memNode->len();

        //first start
        //ST->m_suffixArray[memNode->m_SA_start], length
        //for each start in MEM node
        int *arr = new int[numSeq];
        int flagP = 0, flagNX = 0, pos, l1, l2;
        string mem;
        string leftmem, rightmem;
        vector<int> *memPos1 = new vector<int>[numSeq];
        vector<int> *memPos2 = new vector<int>[numSeq];
        for(int i = 0; i < memNode->m_SA_end - memNode->m_SA_start + 1; i++)
        {
            int y = ST->m_suffixArray[memNode->m_SA_start+i];
            int n = 0, z = 0;

            if(i == 0) {
                mem = S.substr(ST->m_suffixArray[memNode->m_SA_start+i], memNode->m_strdepth);
                if(containsn(mem) || containsx(mem) || containsN(mem) || containsX(mem) ) flagNX = 1;
                //if(containsn(mem) || containsx(mem) || containsN(mem) || containsX(mem) || string_has_all_of_the_same_chars(mem)) flagNX = 1;
                if(containsD(mem)) {
                    flagP = 1;
                    pos = mem.find("#");
                    l1 = pos;
                    l2 = memNode->m_strdepth - pos;
                    leftmem = mem.substr(0, l1);
                    rightmem = mem.substr(pos+1);
                }
            }
            if(flagNX == 1) break;

            if(flagP == 1) {
                if(l1 >= Kmer_Len){
                    seqPos(sLength, y, &n, &z);
                    memPos1[n].push_back(z - 1);
                    //cout<< leftmem << "\t" << "Seq_" << n << "\t" << z-1 << "\t" << ST->m_suffixArray[memNode->m_SA_start+i] << "\t" << l1<<endl;
                }
                if(l2 >= Kmer_Len){
                    seqPos(sLength, y+pos+1, &n, &z);
                    memPos2[n].push_back(z - 1);
                    //cout<< rightmem << "\t" << "Seq_" << n << "\t" << z-1 << "\t" << pos+1 << "\t" << l2<< endl;
                }
            }
            else { // flag == 0
                seqPos(sLength, y, &n, &z);
                memPos1[n].push_back(z - 1);
                //cout << mem << "\t" << "Seq_" << n << "\t" << z-1 << "\t" << ST->m_suffixArray[memNode->m_SA_start+i] << "\t" << memNode->m_strdepth << endl;
            }
        } // end of for loop
        //cout << "testing starts..."<< endl;

        if(flagNX == 0) {
            if(flagP == 0)  {
                cartesian(memPos1, numSeq, memTable, mems, mem);
            }
            if(flagP == 1){
                if(l1 >= Kmer_Len) cartesian(memPos1, numSeq, memTable,mems, leftmem);
                if(l2 >= Kmer_Len) cartesian(memPos2, numSeq, memTable, mems, rightmem);
            }
        }

    }

}

void printMEMs_original(SuffixTree * tree, string S, int *sLength, int numSeq, vector<int *> &memTable, vector<string> &mems)
{
    //recursively perform DFS on suffix tree to find all MEMs
    printMEMnode_original(tree->m_root, tree, S, sLength, numSeq, memTable, mems);

}

void reverseComplemenet(string& s){
//void reverseComplemenet(string& s, int **g, int j){
    /*g[j][0] = s.length() - g[j][0] - 1;
    g[j][1] = s.length() - g[j][1] - 1;*/
    for(int i=0; i<s.length(); i++){
        if(s[i] == 'A') s[i] = 'T';
        else if(s[i] == 'C') s[i] = 'G';
        else if(s[i] == 'T') s[i] = 'A';
        else if(s[i] == 'G') s[i] = 'C';
    }
    //reverse(s, s + s.length());
    reverse(s.begin(), s.end());

}


/* START OF MAIN() Function */

int main__2(int argc, char ** argv)
{
    string S = "";
    bool txt = false;
    bool dot = false;
    bool sort = false;
    bool count = true;
    bool manyMEMs = false;
    bool multiFasta = false;
    bool printModGenome = false;
    int minMEM = 1;

    string filename, kmerLenFile, kvalue;
    string modGenomeFilename;
    string CDGprefix;
    string outf="";

    try{

        for (int c = 1; c < argc; c++)
        {

            if (!strcmp(argv[c], "-h"))       { printHelp(); }
            else if (!strcmp(argv[c], "-mem"))     { MEM = true; kvalue = argv[c+1]; minMEM = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-file"))    { filename = argv[c+1]; c++; }
            else if (!strcmp(argv[c], "-out")) { outf = argv[c+1]; c++; }
        }

        if (filename.empty() || kvalue.empty() || outf.empty())
        {
            printHelp();
        }
        ifstream file;
        file.open(filename.c_str());
        if (! file.is_open())
        {
            cout << "File: " << filename << " does not exist." << endl;
            exit(1);
        }


        SuffixTree * tree;
        cerr << "Loading " << filename << endl;

        //ifstream file;
        //file.open(filename.c_str());
        std::ofstream fastaStartOutfile;

        string buffer;

        if(multiFasta)
        {
            string fileName_str= CDG_Filename + "fastaPos.txt";
            fastaStartOutfile.open(fileName_str.c_str());
        }

        int *sLength;
        int no_seq = 0;

        while(getline(file, buffer))
        {
            //cout<<"length: "<<buffer.length()<<endl;
            //cout<<buffer<<endl<<endl;
            if (buffer[0] == '>')
            {
                if(!isspace(buffer[1])) no_seq++;
                if(!S.empty()) S +='$';


                if(multiFasta)
                {
                    fastaStartOutfile << S.length() <<endl;
                }
            }
            else
            {
                for (int i = 0; i < buffer.length(); i++)
                {
                    char b = toupper(buffer[i]);
                    if (b == 'A' || b == 'C' || b == 'G' || b == 'T' )
                    {
                        S += b;
                    }
                    else //catch ambiguity codes
                    {
                        S += 'N';
                    }
                }
            }
            cerr<<"length:"<<buffer.length()<<"\n";
        }
        //exit(0);
        S += '$';
        sLength = (int *) malloc(sizeof(int) * no_seq);
        printf("The no. of sequences is %d\n", no_seq);
        int max = 0;

        file.close();

        file.open(filename.c_str());
        std::ofstream fastaStartOutfile1;

        if(multiFasta)
        {
            string fileName_str= CDG_Filename + "fastaPos.txt";
            fastaStartOutfile1.open(fileName_str.c_str());
        }


        int i = 0, j = 0, k = 0;

        string *seq = new string[no_seq];
        string *seqName = new string[no_seq];

        int **geneLoc = new int*[no_seq];
        for(int c = 0; c < no_seq; c++)
            geneLoc[c] = new int[2];

        int **newGeneLoc = new int*[no_seq];
        for(int c = 0; c < no_seq; c++)
            newGeneLoc[c] = new int[2];

        string **chrLoc = new string*[no_seq];
        for(int i = 0; i < no_seq; i++)
            chrLoc[i] = new string[4];

        /* For GOBE Visualization files */
        /* ------------------------------
        int x1 = filename.find_last_of("_");
        //int y1 = filename.find_first_of(".");
        int y1 = filename.find_last_of(".");
        string  ID1 = filename.substr(x1+1, y1-x1-1);
        cout<<"ID = " << ID1 <<"\n";
        string  fileGobe= "Gobe"+ID1+".csv";

        char *fileG = new char[fileGobe.length()+1];
        strcpy(fileG, fileGobe.c_str());
        char Gobe[100];
        strcpy(Gobe, fileG);
        delete [] fileG;
        FILE *csv = fopen(Gobe, "w");
        ---------------------------------*/

        string sequence = "";
        while(getline(file, buffer))
        {
            cout<<"\n The length of "<<buffer<<" is :"<<buffer.length()<<"\n";
            if (buffer[0] != '>')
            {
                sequence += buffer;
                //sLength[i]=buffer.length();
                //seq[i]=buffer;
                //i++;
            }
            else if(!isspace(buffer[1]))
            {
                int len = buffer.length();
                //seqName[i] = buffer.substr(1, len);
                //cout<<seqName[i]<<"\n";

                if(i != 0) {
                    seq[i-1] = sequence;
                    sLength[i-1] = seq[i-1].length();
                    sequence = "";
                }
                i++;
            }
            else {
                char gene[buffer.length()+1];
                strcpy(gene, buffer.c_str());
                printf("%s\n", gene);
                int aa = 0;
                char* token = strtok(gene, " \t");
                int start = 0, end = 0;
                while (token) {
                    aa++;
                    printf("token: %s\n", token);
                    if(aa == 2) seqName[k] = token;//geneLoc[k][0] = atoi(token);
                    else if(aa == 3) chrLoc[k][3] = token;
                    else if(aa == 4) chrLoc[k][2] = token;
                    else if(aa == 5) start = atoi(token);
                    else if(aa == 6) chrLoc[k][0] = token;
                    else if(aa == 7) chrLoc[k][1] = token;
                    else if(aa == 8) end = atoi(token);

                    token = strtok(NULL, " \t");
                }
                geneLoc[k][0] = atoi(chrLoc[k][0].c_str()) - start;
                geneLoc[k][1] = atoi(chrLoc[k][1].c_str()) - start;
                cout<<" geneLoc[k][0] "<<  geneLoc[k][0]<< " geneLoc[k][1]"<<geneLoc[k][1]<<endl;
                k++;
            }
        }
        /*cout<<"------------------------\n";
        for(int i=0; i<no_seq; i++)
          for(int j=0; j<2; j++)
            cout<< geneLoc[i][j]<<"\n";
         cout<<"-----------------------\n";  */

        seq[i-1] = sequence;
        sLength[i-1] = seq[i-1].length();
        //exit(0);
        file.close();

        /*string **chrLoc = new string*[no_seq];
        for(int i = 0; i < no_seq; i++)
           chrLoc[i] = new string[4];*/
        //exit(0);
        //getChrLocations(seqName, chrLoc, no_seq, gffname);
        //chrLoc[1][2][0] = '-'; // for testing
        for(int i=0; i<no_seq; i++){
            for(int j=0; j<4; j++)
                cout << chrLoc[i][j]<<" ";
            cout<< "\n";
        }
        //exit(0);
        cout<< "Original gene Locations ..\n";
        for(int c = 0; c < no_seq; c++) {
            newGeneLoc[c][0] = geneLoc[c][0];
            newGeneLoc[c][1] = geneLoc[c][1];
            cout<<newGeneLoc[c][0]<<" "<<newGeneLoc[c][1]<<endl;

        }
        cout<<endl<<endl;
        for(int i=0; i<no_seq; i++){
            if(chrLoc[i][2][0] == '-') {
                //reverseComplemenet(seq[i], geneLoc, i);
                reverseComplemenet(seq[i]);
                cout<<"test:"<<seq[i].length()<<" "<<geneLoc[i][0]<<endl;
                geneLoc[i][0] = seq[i].length() - geneLoc[i][0] - 1;
                geneLoc[i][1] = seq[i].length() - geneLoc[i][1] - 1;
            }
        }

        cout<< " New gene Locations ..\n";
        for(int c = 0; c < no_seq; c++)
            cout<< geneLoc[c][0]<<" "<<geneLoc[c][1]<<"\n";

        i=0;
        /* For GOBE Visualization */
        /*-----------------------------------------------------
        // print the first few lines into Gobe file
        for(j=0; j<no_seq; j++){
           fprintf(csv, "S_%d,S_%d,0,%d,track,", j, j,sLength[j]);
           fprintf(csv,"\n");
        }
        for(i=0; i<no_seq; i++){
          printf("%d %d\n", geneLoc[i][0], geneLoc[i][1]);
          fprintf(csv, "%c%c_Gene,S_%d,%d,%d,gene,+",seqName[i][0],seqName[i][1],i,geneLoc[i][0],geneLoc[i][1]);
          fprintf(csv, "\n");
         }

         fclose(csv);
        -----------------------------------------------------------*/
        j=0;
        for(j = 0; j < no_seq; j++)
            if(max < sLength[j])
                max = sLength[j];

        int total_length = 0;
        for(j = 0; j < no_seq; j++){
            cout<<"length of sequence-"<<j+1<<" is: "<<sLength[j]<<"\n";
            total_length += sLength[j];
            //cout<<"S"<<j+1<<" :"<<seq[j]<<"\n";
        }
        cout << "Total Length "<< total_length << "max ="<<max <<"\n";
        i = 0;
        /* Print the Gene Locations */
        /*----------------------------
        for(i=0; i<no_seq; i++)
          printf("%d %d\n", geneLoc[i][0], geneLoc[i][1]);

        j=0;
        for(j=0; j<i; j++) if(max > sLength[j]) max = sLength[j];
        for(j=0; j<i; j++){
          cout<<"length of "<<j+1<<" sequence is: "<<sLength[j]<<"\n";
          //cout<<"S"<<j+1<<" :"<<seq[j]<<"\n";
        }
        -------------------------------*/
        //exit(0);
        if(multiFasta)
            fastaStartOutfile.close();

        //S += "$";
        S = "s";
        for(i = 0; i < no_seq; i++){
            for(int j = 0; j < seq[i].length(); j++){
                if(!(seq[i][j] == 'A' || seq[i][j] == 'C' || seq[i][j] == 'G' || seq[i][j] == 'T' ))
                    S = S + 'N';
                else S = S + seq[i][j];
            }
            S = S + '#';
        }
        S = S + '$';


        if(printModGenome)
        {
            cerr<<"printing genome string to file"<<endl;
            //save S in file
            ofstream modGenomeFile;
            modGenomeFile.open(modGenomeFilename.c_str());
            modGenomeFile << S;
            modGenomeFile.close();
        }


        cerr << "Creating Suffix Tree for string of length " << S.length() << endl;

        tree = buildUkkonenSuffixTree(S);

        int nodesize = sizeof(SuffixNode);
        treeintLarge totalsize = nodesize * SuffixNode::s_nodecount;

        cerr << "Total nodes: " << SuffixNode::s_nodecount << endl;
        cerr << "Total space: " << totalsize << endl;
        cerr << "Node size: "   << nodesize << endl;

        double bytesperbase = totalsize / ((double) S.length()-2);
        cerr << "Bytes/Base: " << bytesperbase << endl;


        //construct compressed de Bruijn graph for sequence
        int numKmerLens = 0;
        int * kmerLens;

        if(MEM)
        {
            numKmerLens = 1;
            kmerLens = new int[numKmerLens];
            kmerLens[0] = minMEM;
        }
        else if(manyMEMs)
        {
            //store K values in text file, with different K on each line (blank line is ignored)
            // command line argument for manyMEMs that is the filename
            // scan the file twice: once to count number of lines and store in numKmerLens
            // then create array kmerLens = new int[numKmerLens];
            // then scan file a second time to poplate array of kmerLens with actual values of k


            if (!kmerLenFile.empty())
            {
                ifstream kfile;
                string buffer;

                kfile.open(kmerLenFile.c_str());

                while(kfile >> buffer)
                {
                    if(buffer != "")
                        numKmerLens++;
                }
                //return to beginning of file
                kfile.clear();
                kfile.seekg(0, ios::beg);

                kmerLens = new int[numKmerLens];

                int i = 0;
                for(int i = 0; i<numKmerLens; i++)
                {
                    kfile >> buffer;
                    if(buffer != "")
                    {
                        int nextK = atoi(buffer.c_str());
                        kmerLens[i] = nextK;
                    }

                }

                kfile.close();
            }

            CDGprefix = CDG_Filename;
        }

        cerr<<"numKmerLens = "<<numKmerLens<<endl;

        for(int i = 0; i < numKmerLens; i++)
        {
            Kmer_Len = kmerLens[i];
            cerr<<endl<<" Kmer_Len = "<<Kmer_Len<<endl;
            cerr<<" maxMEMstrdepth="<<tree->m_maxMEMstrdepth<<endl;

            if(i == 0)
                tree->markMEMnodes(Kmer_Len, true);
            else
                tree->markMEMnodes(Kmer_Len, false);

            tree->preprocessLMA();
            cerr<<" finished marking MEM nodes"<<endl;

            cerr<<" finished constructing suffix tree and marking MEMs"<<endl;

            vector<int *> memTable;  // vector which contains all starting positions of a MEM
            vector<string> mem_strings; // vector which contains the actual MEM sequence

            //ofstream file1, file2;
            //string fname = "CNS_" + kvalue + ".csv";
            //string mems = "CNSs_" + kvalue + ".txt";
            //file1.open(fname.c_str());
            //file2.open(mems.c_str());
            int cnt = 0;

            printMEMs_original(tree, S, sLength, no_seq, memTable, mem_strings);

            cout<< "size: " << memTable.size() << endl;

            cout<< "size: " << mem_strings.size() << endl;

            int cnt1 = memTable.size();
            int **matrix = new int*[cnt1];
            for(int p = 0; p < cnt1; p++)
                matrix[p] = new int[no_seq];

            int *arr = new int[no_seq];
            string ss;
            i = 0;
            for(std::vector<int *> :: iterator u = memTable.begin(); u != memTable.end(); ++u) {
                arr = *u;
                for(int j = 0; j < no_seq; j++)
                    matrix[i][j] = arr[j];
                i++;
            }


            /*memTable.clear();
            mem_strings.clear();
            
            vector<int *>(memTable).swap(memTable);
            vector<string>(mem_strings).swap(mem_strings);
            */

            int cnt2 = mem_strings.size();
            string  *MEM = new string[cnt2];
            i = 0;
            for(std::vector<string> :: iterator v = mem_strings.begin(); v != mem_strings.end(); ++v) {
                MEM[i] = *v;
                i++;
            }
            cout<<"Testing 1 ..."<<endl;
            delete[] arr;
            int *check = new int[cnt1];
            for(int i=0; i<cnt1; i++) check[i] = 0;

            /* Check if a MEM is subset of other MEM */
            /*for(int p = 0; p < cnt1; p++){
              //int flag1 = 0, flag2 = 0;
              if(check[p] == 0){
                  for(int q = p+1; q <cnt1; q++){
                      int c = 0;
                      if( check[q] == 0){
                          for(int r = 0; r < no_seq; r++){
                              if(matrix[p][r] == matrix[q][r]) c++;
                          }                  
                          if( c == no_seq) {
                              if(MEM[p].length() <= MEM[q].length()) check[p] = 1;
                              else check[q] = 1;
                          } 
                       }  
                    }
                 }
            }*/
            cout<<"Testing 2 ..."<<endl;
            /* Check if CNS is present in all sequences */
            for(int p = 0; p < cnt1; p++){
                int c = 0;
                if(check[p] == 0){
                    for(int q = 0; q < no_seq; q++){
                        if(matrix[p][q] != 0)
                            c++;
                    }
                    if( c != no_seq) check[p] = 1;
                }
            }

            /* Check if the CNS is present on one side of the gene */
            for(int p = 0; p < cnt1; p++){
                int c1 = 0, c2 = 0;
                if(check[p] == 0){
                    for(int q = 0; q < no_seq; q++){
                        if(matrix[p][q] < geneLoc[q][0]) c1++;
                        if(matrix[p][q] > geneLoc[q][1]) c2++;
                    }
                    if(c1 > 0 && c1 < no_seq || c2 > 0 && c2 < no_seq){
                        check[p] = 1;
                    }
                }
            }

            cout<<"Testing 3 ..."<<endl;
            /* Check if the CNS is present in gene region */
            for(int p = 0; p < cnt1; p++){
                int c1 = 0, c2 = 0;
                if(check[p] == 0){
                    for(int q = 0; q < no_seq; q++){
                        if(matrix[p][q] >= geneLoc[q][0] && matrix[p][q] <= geneLoc[q][1]) c1++;
                        if(matrix[p][q] <= geneLoc[q][0] && matrix[p][q] >= geneLoc[q][1]) c2++;
                    }
                    if(c1 != 0) check[p] = 1;
                    if(c2 != 0) check[p] = 1;
                }
            }


            /* moved to here ... on 5/21/2016 */
            /*for(int p = 0; p < cnt1; p++){
              //int flag1 = 0, flag2 = 0;
              if(check[p] == 0){
                  for(int q = p+1; q <cnt1; q++){
                      int c = 0;
                      if( check[q] == 0){
                          for(int r = 0; r < no_seq; r++){
                              if(matrix[p][r] == matrix[q][r]) c++;
                          }                  
                          if( c == no_seq) {
                              if(MEM[p].length() <= MEM[q].length()) check[p] = 1;
                              else check[q] = 1;
                          } 
                       }  
                    }
                 }
            }
            
            */

            int size = 0;
            for(int i=0; i<cnt1; i++)
                if(check[i] == 0) size++;
            cout<< "size= "<<size<<endl;
            int **new_matrix = new int*[size];
            for(int i=0; i<size; i++)
                new_matrix[i] =  new int[no_seq];

            string *new_mems = new string[size];

            int ind=0;
            for(int i=0; i<cnt1; i++){
                if(check[i] == 0) {
                    new_mems[ind] = MEM[i];
                    for(int j=0; j<no_seq; j++)
                        new_matrix[ind][j] = matrix[i][j];
                    ind++;
                }
            }
            cout<<"Testing 4 ..."<<endl;
            /* moved to here ... on 5/21/2016 */
            for(int p = 0; p < size; p++){
                if(p%10000 == 0) cout<< "10000 completed.."<<endl;
                for(int q = p+1; q <size; q++){
                    int c = 0;
                    for(int r = 0; r < no_seq; r++){
                        if(new_matrix[p][r] == new_matrix[q][r]) c++;
                    }
                    if( c == no_seq) {
                        if(new_mems[p].length() <= new_mems[q].length()) check[p] = 1;
                        else check[q] = 1;
                    }
                }
            }

            cout<<"Testing 5 ..."<<endl;
            //exit(0);
            ofstream file1, file2, file3, file4, file5;
            outf += "_";
            string  fileM1= outf + "mems_" + kvalue + ".txt";
            string  fileM2= outf + "MEM_1_" +kvalue + ".csv";
            string  fileM3= outf + "MEM_2_" +kvalue + ".csv";

            file1.open(fileM1.c_str());
            file2.open(fileM2.c_str());
            file3.open(fileM3.c_str());

            //cout<<"______________________________________"<<"\n";
            file2<< "MEM,Length";
            file3<< "Length";
            for(int i=0; i<no_seq; i++){
                file2<< ","<<"S_"<< (i+1);
                file3<< ","<<"S_"<< (i+1);
            }
            file2<<"\n";
            file3<<"\n";

            //cout<<"ID\tLength\t";
            for(int p = 0; p < no_seq; p++)
                //cout<<"S_"<<p+1<<"\t";
                //cout<<"\n";
                file2 << "0";
            file3 << "0";
            for(int p = 0; p < no_seq; p++){
                file2 << ",0";
                file3 << ",0";
            }
            file2<<"\n";
            file3<<"\n";

            int l = 1;
            for(int p = 0; p < cnt1; p++){
                if(check[p] == 0 ){
                    //cout << l;
                    file2 << l;
                    //file3 << "," << l;
                    //cout << "\t"  << MEM[p].length();
                    file2 << "," << MEM[p].length();
                    file3 << MEM[p].length();
                    file1 << MEM[p] << "\n";

                    for(int q = 0; q < no_seq; q++){
                        // cout << "\t" << matrix[p][q];
                        file2 << "," << matrix[p][q];
                        file3 << "," << matrix[p][q];
                    }
                    l++;
                    //cout <<"\n";
                    file2<<"\n";
                    file3<<"\n";
                }
            }
            //}
            //cout<< "Max = "<< max<<endl;
            file2 << "0";
            file3 << "0";
            for(int p = 0; p < no_seq; p++){
                file2 << ","<<max;
                file3 << ","<<max;
            }


            //cout<<"testing again2\n";
            //exit(0);
            cout<<"______________________________________"<<"\n";
            file1.close();
            file2.close();
            file3.close();
            delete[] matrix;
            delete[] MEM;

            cout<< seqName[0]<<endl;
            vector<int *> LP_Matrix;
            vector<string> LP_MEM;
            cout<< "Calling LP .."<<endl;
            // int **new_matrix = new int*[size];
            // string *new_mems = new string[size];
            // new_mems[ind] = MEM[i];
            // string *seqName = new string[no_seq];
            LP(new_matrix, new_mems, size, no_seq, LP_Matrix, LP_MEM, seqName);   // here is the function that I want a lot
            cout<<"after calling LP()..."<<endl;
            //delete [] sLength;
            //delete [] seqName;
            //delete [] geneLoc;
            //delete [] MEM;
            //delete [] matrix;
            //cout<< seqName[0]<<endl;
            //exit(0);

            string fileM4 = outf + "CNS_" +kvalue + ".csv";
            string fileM5 = outf + "LPMEM_" + kvalue + ".csv";

            file5.open(fileM5.c_str());
            file5<< "Length";
            for(int i=0; i<no_seq; i++){
                file5<< ","<<"S_"<< (i+1);
            }

            file5<<"\n";
            file5 << "0";
            for(int p = 0; p < no_seq; p++){
                file5 << ",0";
            }
            file5<<"\n";


            file4.open(fileM4.c_str());
            //	cout<< "Testing before Exit: "<< seqName[i].c_str() <<endl;
            //exit(0);
            for(int i=0; i<no_seq; i++)
                file4 << seqName[i].c_str()<<"\n";

            file4 << "LEN";


            for(int i=0; i<no_seq; i++)
                file4 << ","<<seqName[i][0] << seqName[i][1];
            for(int i=0; i<no_seq; i++)
                file4 << ",Chr_"<<seqName[i][0] << seqName[i][1];
            file4 << ",CNS";
            file4 << "\n";

            int *arr1 = new int[no_seq];
            i = 0;
            for(std::vector<int *> :: iterator u = LP_Matrix.begin(); u != LP_Matrix.end(); ++u) {
                arr1 = *u;
                int x2 = LP_MEM.at(i).length();
                file5 << x2;
                file4 << x2;
                //i++;
                for(int j = 0; j < no_seq; j++){
                    int x3 = arr1[j];
                    file5<< ","<<x3;
                    //if(chrLoc[0][no_seq-1][0] == '+'){
                    if(chrLoc[j][2][0] == '+'){
                        x3 = x3 + atoi(chrLoc[j][0].c_str()) - newGeneLoc[j][0];
                    }
                    else{
                        x3 = atoi(chrLoc[j][1].c_str()) + (seq[j].length() - newGeneLoc[j][1]) - x3 - x2;
                    }
                    file4<< ","<<x3;
                }
                for(int j = 0; j < no_seq; j++)
                    file4<< ","<< chrLoc[j][3].c_str();
                file4 << ","<<LP_MEM.at(i);
                file4 << "\n";
                file5<< "\n";
                i++;

            }
            file5 << "0";
            for(int p = 0; p < no_seq; p++){
                file5 << ","<<max;
            }
            file4.close();
            file5.close();

            delete[] new_matrix;
            delete[] new_mems;
            delete[] seq;
            delete[] geneLoc;
            delete[] newGeneLoc;

            /* For HTML Visualization file */
            ifstream vizIn("Sample.html");
            ofstream vizOut1, vizOut2;
            string  html1 = outf + "MEM_"+ kvalue + ".html";
            string  html2 = outf + "LPMEM_"+ kvalue + ".html";
            vizOut1.open(html1.c_str());
            vizOut2.open(html2.c_str());
            /*
            var spec = ["S_1", "S_2","S_3", "S_4",	"S_5"]; // CHANGE-1
            var XC = [0, 0, 0, 0, 0]; // CHANGE-2
            var YC = [0, 0, 0, 0, 0]; // CHANGE-3 
            var n = 5;// CHANGE-4
            d3.csv("MEM_30.csv", function(error, cars) {
            */
            string line;
            if (vizIn.is_open())
            {
                while ( getline (vizIn, line) )
                {
                    if(line.find("CHANGE-1") != string::npos){
                        vizOut1 << "var spec = [\"S_1\"";
                        vizOut2 << "var spec = [\"S_1\"";
                        for(int i=1; i<no_seq; i++){
                            vizOut1 << ", \"S_"<<(i+1)<<"\"";
                            vizOut2 << ", \"S_"<<(i+1)<<"\"";
                        }
                        vizOut1 << "];\n";
                        vizOut2 << "];\n";
                    }
                    else if(line.find("CHANGE-2") != string::npos){
                        vizOut1 << "var XC = [ 0";
                        vizOut2 << "var XC = [ 0";
                        for(int i=0; i<no_seq-1; i++){
                            vizOut1 << ", 0";
                            vizOut2 << ", 0";
                        }
                        vizOut1 << "];\n";
                        vizOut2 << "];\n";
                    }
                    else if(line.find("CHANGE-3") != string::npos){
                        vizOut1 << "var YC = [ 0";
                        vizOut2 << "var YC = [ 0";
                        for(int i=0; i<no_seq-1; i++){
                            vizOut1 << ", 0";
                            vizOut2 << ", 0";
                        }
                        vizOut1 << "];\n";
                        vizOut2 << "];\n";
                    }
                    else if(line.find("CHANGE-4") != string::npos){
                        vizOut1 << "var n = "<< no_seq <<";\n";
                        vizOut2 << "var n = "<< no_seq <<";\n";
                    }
                    else if(line.find("CHANGE-5") != string::npos){
                        vizOut1 << "d3.csv(\""<<fileM3<<"\", function(error, cars) {"<<"\n";
                        vizOut2 << "d3.csv(\""<<fileM5<<"\", function(error, cars) {"<<"\n";
                    }
                    else{
                        vizOut1 << line << '\n';
                        vizOut2 << line << '\n';
                    }
                }
                vizIn.close();
                vizOut1.close();
                vizOut2.close();
            }



        }

        delete [] kmerLens;
        /*delete[] seq;
        delete[] geneLoc;
        delete[] newGeneLoc;*/


        if (txt)  { tree->dumpTreeText(cout); }
        else if (dot)  { tree->dumpTree(); cerr<<" dumping ST to dot file";}
        else if (sort) { tree->dumpTreeSorted(cout, tree->m_root, ""); }


    }
    catch(std::bad_alloc& ba){
        cerr<<"bad_alloc caught "<<ba.what()<<endl;

    }
    catch(exception& ex){
        cerr<<" some other exception caught"<<ex.what() <<endl;
    }

}
