//#include <mmintrin.h>
//#include <xmmintrin.h> 

#include <vector>
#include <cmath>
//#include <omp.h>
//#include <algorithm>
//#include <execution>
#include <iostream>
#include <vector>


#include "algorithm"
#include "execution"

namespace TSnap {


/////////////////////////////////////////////////
// Node centrality measures
double GetDegreeCentr(const PUNGraph& Graph, const int& NId) {
  if (Graph->GetNodes() > 1) {
    return double(Graph->GetNI(NId).GetDeg())/double(Graph->GetNodes()-1); }
  else { return 0.0; }
}

void GetEigenVectorCentr(const PUNGraph& Graph, TIntFltH& NIdEigenH, const double& Eps, const int& MaxIter) {
  const int NNodes = Graph->GetNodes();
  NIdEigenH.Gen(NNodes);
  // initialize vector values
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    NIdEigenH.AddDat(NI.GetId(), 1.0/NNodes);
    IAssert(NI.GetId() == NIdEigenH.GetKey(NIdEigenH.Len()-1));
  }
  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {
    int j = 0;
    // add neighbor values
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      TmpV[j] = 0;
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        TmpV[j] += NIdEigenH.GetDat(NI.GetOutNId(e)); }
    }

    // normalize
    double sum = 0;
    for (int i = 0; i < TmpV.Len(); i++) {
      sum += (TmpV[i]*TmpV[i]);
    }
    sum = sqrt(sum);
    for (int i = 0; i < TmpV.Len(); i++) {
      TmpV[i] /= sum;
    }

    // compute difference
    double diff = 0.0;
    j = 0;
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      diff += fabs(NIdEigenH.GetDat(NI.GetId())-TmpV[j]);
    }

    // set new values
    j = 0;
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      NIdEigenH.AddDat(NI.GetId(), TmpV[j]);
    }

    if (diff < Eps) {
      break;
    }
  }
}

// Group centrality measures
double GetGroupDegreeCentr(const PUNGraph& Graph, const PUNGraph& Group) {
  int deg;
  TIntH NN;
  for (TUNGraph::TNodeI NI = Group->BegNI(); NI < Group->EndNI(); NI++) { 
    deg = Graph->GetNI(NI.GetId()).GetDeg();
    for (int i=0; i<deg; i++) {
      if (Group->IsNode(Graph->GetNI(NI.GetId()).GetNbrNId(i))==0)
      NN.AddDat(Graph->GetNI(NI.GetId()).GetNbrNId(i),NI.GetId());
    }
  }
  return (double)NN.Len();
}

double GetGroupDegreeCentr0(const PUNGraph& Graph, const TIntH& GroupNodes) {
  int deg;
  TIntH NN;
  for (int i = 0; i<GroupNodes.Len(); i++) { 
    deg = Graph->GetNI(GroupNodes.GetDat(i)).GetDeg();
    for (int j = 0; j < deg; j++) {
      if (GroupNodes.IsKey(Graph->GetNI(GroupNodes.GetDat(i)).GetNbrNId(j))==0)
      NN.AddDat(Graph->GetNI(GroupNodes.GetDat(i)).GetNbrNId(j),GroupNodes.GetDat(i));
    }
  }
  return (double)NN.Len();
}

double GetGroupDegreeCentr(const PUNGraph& Graph, const TIntH& GroupNodes) {
  int deg;
  TIntH NN;
  TIntH GroupNodes1;

  for (THashKeyDatI<TInt,TInt> NI = GroupNodes.BegI(); NI < GroupNodes.EndI(); NI++)
    GroupNodes1.AddDat(NI.GetDat(),NI.GetDat());

  for (THashKeyDatI<TInt,TInt> NI = GroupNodes1.BegI(); NI < GroupNodes1.EndI(); NI++){
    TUNGraph::TNodeI node = Graph->GetNI(NI.GetKey());
    deg = node.GetDeg();
    for (int j = 0; j < deg; j++){
      if (GroupNodes1.IsKey(node.GetNbrNId(j))==0 && NN.IsKey(node.GetNbrNId(j))==0)
        NN.AddDat(node.GetNbrNId(j),NI.GetKey());
    }
  }

  return (double)NN.Len();
}

double GetGroupFarnessCentr(const PUNGraph& Graph, const TIntH& GroupNodes) {
  TIntH* NDistH = new TIntH[GroupNodes.Len()];
  
  for (int i=0; i<GroupNodes.Len(); i++){
    NDistH[i](Graph->GetNodes());
  TSnap::GetShortPath<PUNGraph>(Graph, GroupNodes.GetDat(i), NDistH[i], true, TInt::Mx);
  }

  int min, dist, sum=0, len=0;
  for (PUNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    if(NDistH[0].IsKey(NI.GetId()))
      min = NDistH[0].GetDat(NI.GetId());
    else
      min = -1;
    for (int j=1; j<GroupNodes.Len(); j++){
    if (NDistH[j].IsKey(NI.GetId()))
      dist = NDistH[j].GetDat(NI.GetId());
      else
      dist = -1;
    if ((dist < min && dist != -1) || (dist > min && min == -1))
      min = dist;
    }
    if (min>0){  
      sum += min;
      len++;
  }
    
  }

  if (len > 0) { return sum/double(len); }
  else { return 0.0; }
}

PUNGraph *AllGraphsWithNNodes(int n){
  PUNGraph* g = new PUNGraph[(((n*n)-n)/2)+1];
  PUNGraph g0;
  for(int i=0; i<n; i++)
    g0->AddNode(i);
  
  g[0] = g0;
  int br=1;

  for(int i=0; i<n; i++)
  for(int j=i; j<n; j++){
      g0->AddEdge(i,j);
    g[br] = g0;
    br++;
  }

  return g;
}

TIntH *AllCombinationsMN(int m, int n){
  float N = 1;
  for(int i=n; i>0; i--){
    N *= (float)m/(float)n;
    m--;
  n--;
  }

  TIntH* C = new TIntH[(int)N];
  return C;
}

double GetGroupClosenessCentr(const PUNGraph& Graph, const TIntH& GroupNodes) {
  const double Farness = GetGroupFarnessCentr(Graph, GroupNodes);
  if (Farness != 0.0) { return 1.0/Farness; }
  else { return 0.0; }
}

TIntH MaxCPGreedyBetter(const PUNGraph& Graph, const int k) {
  TIntH GroupNodes; // buildup cpntainer of group nodes
  TIntH NNodes; // container of neighbouring nodes
  TIntH Nodes; // nodes sorted by vd
  double gc = 0, gc0 = 0;
  int addId = 0, addIdPrev = 0;
  
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    Nodes.AddDat(NI.GetId(),NI.GetDeg());
  }

  Nodes.SortByDat(false);

  int br = 0;
  while (br < k) {
    for (THashKeyDatI<TInt,TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++) {
      if ((NI.GetDat() <= (int)gc0))
        break;
      gc = NI.GetDat()-Intersect(Graph->GetNI(NI.GetKey()),NNodes);
      if (gc>gc0) {
        gc0 = gc;
        addId = NI.GetKey();
      }
    }
  
    if (addId != addIdPrev){

      GroupNodes.AddDat(br,addId);
      br++;
      gc0=0;

      NNodes.AddDat(addId,0);
      for (int i=0; i<Graph->GetNI(addId).GetDeg(); i++) {
        NNodes.AddDat(Graph->GetNI(addId).GetNbrNId(i),0);
      }
      addIdPrev = addId;
      Nodes.DelKey(addId);
    } else {
      br = k;
    }
    printf("%i,",br);
  }

  // gcFinal = GetGroupDegreeCentr(Graph, GroupNodes);
  return GroupNodes;
}

// this is the variation of the first version that doesent stop after finding the optimal K
TIntH MaxCPGreedyBetter1(const PUNGraph& Graph, const int k) {
  TIntH GroupNodes;
  TIntH NNodes;
  TIntH Nodes;
  double gc = 0, gc0 = 0;
  int addId = 0, addIdPrev = 0;
  
  // put nodes in the container and sort them by vertex degree
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    Nodes.AddDat(NI.GetId(),NI.GetDeg());
  }
  Nodes.SortByDat(false);

  int br = 0;
  while (br < k) {
    for (THashKeyDatI<TInt,TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++){
      if((NI.GetDat() < (int)gc0))
        break;
      gc = NI.GetDat()-Intersect(Graph->GetNI(NI.GetKey()),NNodes);
      if (gc>gc0) {
        gc0 = gc;
        addId = NI.GetKey();
      }
    }
  
    if (addId != addIdPrev){

      GroupNodes.AddDat(br,addId);
      br++;
      gc0=-10000000;
  
      NNodes.AddDat(addId,0);
      for (int i=0; i<Graph->GetNI(addId).GetDeg(); i++) {
        NNodes.AddDat(Graph->GetNI(addId).GetNbrNId(i),0);
      }
      addIdPrev = addId;
      Nodes.DelKey(addId);
    }
  }

  // gcFinal = GetGroupDegreeCentr(Graph, GroupNodes);
  return GroupNodes;
}

// version with string type of container of group nodes - Fail (it is slower)
TIntH MaxCPGreedyBetter2(const PUNGraph& Graph, const int k) {
  TIntH GroupNodes; // buildup cpntainer of group nodes
  TStr NNodes; // container of neighbouring nodes
  TIntH Nodes; // nodes sorted by vd
  double gc = 0, gc0 = 0;
  int addId = 0, addIdPrev=0;
  
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    Nodes.AddDat(NI.GetId(),NI.GetDeg());
  }

  Nodes.SortByDat(false);

  int br=0;
  while (br < k) {
    for (THashKeyDatI<TInt,TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++){
      if((NI.GetDat() <= (int)gc0))
        break;
      gc = NI.GetDat()-Intersect(Graph->GetNI(NI.GetKey()),NNodes);
      if (gc>gc0) {
        gc0 = gc;
        addId = NI.GetKey();
      }
    }
  
    if (addId != addIdPrev) {

      GroupNodes.AddDat(br,addId);
      br++;
      gc0=0;
    
      TInt digi = addId;
      TStr buf = digi.GetStr();

      NNodes += " "+buf;

      for (int i=0; i<Graph->GetNI(addId).GetDeg(); i++) {
        TInt digi = Graph->GetNI(addId).GetNbrNId(i);
        TStr buf = digi.GetStr();
        NNodes += " "+buf;
      }
      addIdPrev = addId;
      Nodes.DelKey(addId);
    } else {
      br = k;
    }
    printf("%i,",br);
  }

  // gcFinal = GetGroupDegreeCentr(Graph, GroupNodes);
  return GroupNodes;
}

// version with int array - the fastest
TIntH MaxCPGreedyBetter3(const PUNGraph& Graph, const int k) {
  TIntH GroupNodes; // buildup cpntainer of group nodes
  const int n = Graph->GetNodes();
  int *NNodes = new int[n]; // container of neighbouring nodes
  int NNodes_br = 0;
  TIntH Nodes; // nodes sorted by vd
  double gc = 0, gc0 = 0;
  int addId = 0, addIdPrev = 0;
  
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    Nodes.AddDat(NI.GetId(),NI.GetDeg());
  }

  Nodes.SortByDat(false);

  int br = 0;
  while (br < k) {
    for (THashKeyDatI<TInt,TInt> NI = Nodes.BegI(); NI < Nodes.EndI(); NI++){
      if((NI.GetDat() <= (int)gc0))
        break;
      gc = NI.GetDat()-Intersect(Graph->GetNI(NI.GetKey()),NNodes,NNodes_br);
      if (gc>gc0){
        gc0 = gc;
        addId = NI.GetKey();
      }
    }
  
    if (addId != addIdPrev) {

      GroupNodes.AddDat(br,addId);
      br++;
      gc0=0;

      int nn = addId;
      bool nnnew = true;
      for (int j=0; j<NNodes_br; j++)
        if (NNodes[j] == nn){
          nnnew = false;
          j = NNodes_br;
        }
  
      if (nnnew){
        NNodes[NNodes_br] = nn;
        NNodes_br++;
      }

      for (int i=0; i<Graph->GetNI(addId).GetDeg(); i++) {
        int nn = Graph->GetNI(addId).GetNbrNId(i);
        bool nnnew = true;
        for (int j=0; j<NNodes_br; j++) {
          if (NNodes[j] == nn){
            nnnew = false;
            j = NNodes_br;
          }
        }
        if (nnnew){
          NNodes[NNodes_br] = nn;
          NNodes_br++;
        }
      }
      addIdPrev = addId;
      Nodes.DelKey(addId);
    } else {
      br = k;
    }
    printf("%i,",br);
  }

  delete NNodes;
  // gcFinal = GetGroupDegreeCentr(Graph, GroupNodes);
  return GroupNodes;
}

//Weighted PageRank
int GetWeightedPageRank(const PNEANet Graph, TIntFltH& PRankH, const TStr& Attr, const double& C, const double& Eps, const int& MaxIter) {
  if (!Graph->IsFltAttrE(Attr)) return -1;

  TFltV Weights = Graph->GetFltAttrVecE(Attr);

  int mxid = Graph->GetMxNId();
  TFltV OutWeights(mxid);
  Graph->GetWeightOutEdgesV(OutWeights, Weights);

  const int NNodes = Graph->GetNodes();
  //const double OneOver = 1.0/double(NNodes);
  PRankH.Gen(NNodes);
  for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    PRankH.AddDat(NI.GetId(), 1.0/NNodes);
    //IAssert(NI.GetId() == PRankH.GetKey(PRankH.Len()-1));
  }
  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {
    int j = 0;
    for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
      TmpV[j] = 0;
      for (int e = 0; e < NI.GetInDeg(); e++) {
        const int InNId = NI.GetInNId(e);
        const TFlt OutWeight = OutWeights[InNId];
        int EId = Graph->GetEId(InNId, NI.GetId());
        const TFlt Weight = Weights[Graph->GetFltKeyIdE(EId)];
        if (OutWeight > 0) {
          TmpV[j] += PRankH.GetDat(InNId) * Weight / OutWeight; }
      }
      TmpV[j] =  C*TmpV[j]; // Berkhin (the correct way of doing it)
      //TmpV[j] =  C*TmpV[j] + (1.0-C)*OneOver; // iGraph
    }
    double diff=0, sum=0, NewVal;
    for (int i = 0; i < TmpV.Len(); i++) { sum += TmpV[i]; }
    const double Leaked = (1.0-sum) / double(NNodes);
    for (int i = 0; i < PRankH.Len(); i++) { // re-instert leaked PageRank
      NewVal = TmpV[i] + Leaked; // Berkhin
      //NewVal = TmpV[i] / sum;  // iGraph
      diff += fabs(NewVal-PRankH[i]);
      PRankH[i] = NewVal;
    }
    if (diff < Eps) { break; }
  }
  return 0;
}

#ifdef USE_OPENMP
int GetWeightedPageRankMP(const PNEANet Graph, TIntFltH& PRankH, const TStr& Attr, const double& C, const double& Eps, const int& MaxIter) {
  if (!Graph->IsFltAttrE(Attr)) return -1;
  const int NNodes = Graph->GetNodes();
  TVec<TNEANet::TNodeI> NV;

  //const double OneOver = 1.0/double(NNodes);
  PRankH.Gen(NNodes);
  int MxId = 0;

  for (TNEANet::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    NV.Add(NI);
    PRankH.AddDat(NI.GetId(), 1.0/NNodes);
    int Id = NI.GetId();
    if (Id > MxId) {
      MxId = Id;
    }
  }

  TFltV PRankV(MxId+1);
  TFltV OutWeights(MxId+1);

  TFltV Weights = Graph->GetFltAttrVecE(Attr);

  #pragma omp parallel for schedule(dynamic,10000)
  for (int j = 0; j < NNodes; j++) {
    TNEANet::TNodeI NI = NV[j];
    int Id = NI.GetId();
    OutWeights[Id] = Graph->GetWeightOutEdges(NI, Attr);
    PRankV[Id] = 1/NNodes;
  }

  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {

    #pragma omp parallel for schedule(dynamic,10000)
    for (int j = 0; j < NNodes; j++) {
      TNEANet::TNodeI NI = NV[j];
      TFlt Tmp = 0;
      for (int e = 0; e < NI.GetInDeg(); e++) {
        const int InNId = NI.GetInNId(e);

        const TFlt OutWeight = OutWeights[InNId];

        int EId = Graph->GetEId(InNId, NI.GetId());
        const TFlt Weight = Weights[Graph->GetFltKeyIdE(EId)];

        if (OutWeight > 0) {
          Tmp += PRankH.GetDat(InNId) * Weight / OutWeight;
        }
      }
      TmpV[j] =  C*Tmp; // Berkhin (the correct way of doing it)
      //TmpV[j] =  C*TmpV[j] + (1.0-C)*OneOver; // iGraph
    }

    double sum = 0;
    #pragma omp parallel for reduction(+:sum) schedule(dynamic,10000)
    for (int i = 0; i < TmpV.Len(); i++) { sum += TmpV[i]; }
    const double Leaked = (1.0-sum) / double(NNodes);

    double diff = 0;
    #pragma omp parallel for reduction(+:diff) schedule(dynamic,10000)
    for (int i = 0; i < NNodes; i++) {
      TNEANet::TNodeI NI = NV[i];
      double NewVal = TmpV[i] + Leaked; // Berkhin
      //NewVal = TmpV[i] / sum;  // iGraph
      int Id = NI.GetId();
      diff += fabs(NewVal-PRankV[Id]);
      PRankV[Id] = NewVal;
    }
    if (diff < Eps) { break; }
  }

  #pragma omp parallel for schedule(dynamic,10000)
  for (int i = 0; i < NNodes; i++) {
    TNEANet::TNodeI NI = NV[i];
    PRankH[i] = PRankV[NI.GetId()];
  }

  return 0;
}

#endif // USE_OPENMP

//Event importance
TIntFltH EventImportance(const PNGraph& Graph, const int k) {
  TIntFltH NodeList; // values for nodese

  for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    NodeList.AddDat(NI.GetId(),NI.GetOutDeg());
  }


  for (THashKeyDatI<TInt,TFlt> NI = NodeList.BegI(); NI < NodeList.EndI(); NI++){
    int outdeg = Graph->GetNI(NI.GetKey()).GetOutDeg();
    int indeg = Graph->GetNI(NI.GetKey()).GetInDeg();
    
    if (outdeg>1 && indeg>0){
      double val = (1-(1/(double)outdeg))/(double)indeg;
      for(int i=0; i<(outdeg+indeg);i++){
        int NId = Graph->GetNI(NI.GetKey()).GetNbrNId(i);
        if (Graph->GetNI(NI.GetKey()).IsInNId(NId) == true){
        NodeList.AddDat(NId,NodeList.GetDat(NId)+val);
        }
        
      }
    }
    
  }

  return NodeList;
}

//Event importance 1
TIntFltH EventImportance1 (const PNGraph& Graph, const int k) {
  TIntFltH NodeList; // values for nodese

  for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++){
    NodeList.AddDat(NI.GetId(),NI.GetOutDeg());
  }


  for (THashKeyDatI<TInt,TFlt> NI = NodeList.BegI(); NI < NodeList.EndI(); NI++){
    int outdeg = Graph->GetNI(NI.GetKey()).GetOutDeg();
    int indeg = Graph->GetNI(NI.GetKey()).GetInDeg();
    
    if (outdeg>1 && indeg>0){
      double val = (1-(1/(double)outdeg))/(double)indeg;
      for(int i=0; i<(outdeg+indeg);i++){
        int NId = Graph->GetNI(NI.GetKey()).GetNbrNId(i);
        if (Graph->GetNI(NI.GetKey()).IsInNId(NId) == true){
        NodeList.AddDat(NId,NodeList.GetDat(NId)+val);
        }
        
      }
    }
    
  }

  return NodeList;
}

int Intersect(TUNGraph::TNodeI Node, TIntH NNodes){
  int br=0;
  for (int i=0; i<Node.GetDeg(); i++)
  {
    if (NNodes.IsKey(Node.GetNbrNId(i)))
      br++;
  }
  if (NNodes.IsKey(Node.GetId()))
    br++;

  return br;
}

int Intersect(TUNGraph::TNodeI Node, TStr NNodes){
  int br=0;

  TInt digi = -1;
  TStr buf = "";

  for (int i=0; i<Node.GetDeg(); i++)
  {
    digi = Node.GetNbrNId(i);
    TStr buf = digi.GetStr();

    if (NNodes.IsStrIn(buf.CStr()))
    br++;
  }

  digi = Node.GetId();
  buf = digi.GetStr();

  if (NNodes.IsStrIn(buf.CStr()))
    br++;

  return br;
}

int Intersect(TUNGraph::TNodeI Node, int *NNodes, int NNodes_br){
  int br = 0;
  int neig;
  for (int i=0; i<Node.GetDeg(); i++)
  {
    neig = Node.GetNbrNId(i);
    for (int j=0; j<NNodes_br; j++)
    {
    if (neig == NNodes[j])
    {
      br++;
      j = NNodes_br;
    }
    }
  }

  neig = Node.GetId();
  for (int j=0; j<NNodes_br; j++)
  {
    if (neig == NNodes[j])
    {
      br++;
      j = NNodes_br;
    }
  }

  return br;
}

int Intersect1(TUNGraph::TNodeI Node, TStr NNodes){
  int br=0;
  for (int i=0; i<Node.GetDeg(); i++)
  {
    TInt digi = Node.GetNbrNId(i);
    TStr buf = "";
    buf = digi.GetStr();

    if (NNodes.SearchStr(buf.CStr())!=-1)
    br++;
  }
  
  TInt digi = Node.GetId();
  TStr buf = digi.GetStr();

  if (NNodes.SearchStr(buf.CStr())!=-1)
    br++;

  return br;
}

TIntH LoadNodeList(TStr InFNmNodes){
  TSsParser Ss(InFNmNodes, ssfWhiteSep, true, true, true);
  TIntIntH Nodes;
  int br = 0, NId;
  while (Ss.Next()) {
    if (Ss.GetInt(0, NId)) { 
    Nodes.AddDat(br,NId);
    br++;
  }
  }
  return Nodes;
}


int findMinimum(TIntV& Frontier, TIntFltH& NIdDistH) {
  TFlt minimum = TInt::Mx;
  int min_index = 0;
  for (int i = 0; i < Frontier.Len(); i++) {
    int NId = Frontier.GetVal(i);
    if (NIdDistH.GetDat(NId) < minimum) {
      minimum = NIdDistH.GetDat(NId);
      min_index = i;
    }
  }
  const int NId = Frontier.GetVal(min_index);
  Frontier.Del(min_index);
  return NId;
}

int GetWeightedShortestPath(
const PNEANet Graph, const int& SrcNId, TIntFltH& NIdDistH, const TFltV& Attr) {
  TIntV frontier;

  NIdDistH.Clr(false); NIdDistH.AddDat(SrcNId, 0);
  frontier.Add(SrcNId);
  while (! frontier.Empty()) {
    const int NId = findMinimum(frontier, NIdDistH);
    const PNEANet::TObj::TNodeI NodeI = Graph->GetNI(NId);
    for (int v = 0; v < NodeI.GetOutDeg(); v++) {
      int DstNId = NodeI.GetOutNId(v);
      int EId = NodeI.GetOutEId(v);

      if (! NIdDistH.IsKey(DstNId)) {
        NIdDistH.AddDat(DstNId, NIdDistH.GetDat(NId) + Attr[EId]);
        frontier.Add(DstNId);
      } else {
        if (NIdDistH.GetDat(DstNId) > NIdDistH.GetDat(NId) + Attr[EId]) {
          NIdDistH.GetDat(DstNId) = NIdDistH.GetDat(NId) + Attr[EId]; 
        }
      }
    }
  }
  return 0;
}

double GetWeightedFarnessCentr(const PNEANet Graph, const int& NId, const TFltV& Attr, const bool& Normalized, const bool& IsDir) {
  TIntFltH NDistH(Graph->GetNodes());
  
  GetWeightedShortestPath(Graph, NId, NDistH, Attr);
  
  double sum = 0;
  for (TIntFltH::TIter I = NDistH.BegI(); I < NDistH.EndI(); I++) {
    sum += I->Dat();
  }
  if (NDistH.Len() > 1) { 
    double centr = sum/double(NDistH.Len()-1); 
    if (Normalized) {
      centr *= (Graph->GetNodes() - 1)/double(NDistH.Len()-1);
    }
    return centr;
  }
  else { return 0.0; }
}

double GetWeightedClosenessCentr(const PNEANet Graph, const int& NId, const TFltV& Attr, const bool& Normalized, const bool& IsDir) {
  const double Farness = GetWeightedFarnessCentr(Graph, NId, Attr, Normalized, IsDir);
  if (Farness != 0.0) { return 1.0/Farness; }
  else { return 0.0; }
  return 0.0;
}

void GetWeightedBetweennessCentr(const PNEANet Graph, const TIntV& BtwNIdV, TIntFltH& NodeBtwH, const bool& DoNodeCent, TIntPrFltH& EdgeBtwH, const bool& DoEdgeCent, const TFltV& Attr, const bool& IsDir) {
  if (DoNodeCent) { NodeBtwH.Clr(); }
  if (DoEdgeCent) { EdgeBtwH.Clr(); }
  const int nodes = Graph->GetNodes();
  TIntS S(nodes);
  TIntQ Q(nodes);
  TIntIntVH P(nodes); // one vector for every node
  TIntFltH delta(nodes);
  TIntFltH sigma(nodes), d(nodes);
  // init
  for (PNEANet::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    if (DoNodeCent) {
      NodeBtwH.AddDat(NI.GetId(), 0); }
    if (DoEdgeCent) {
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        if (Graph->HasFlag(gfDirected) && IsDir) {
          // add all outgoing edges for directed graphs
          EdgeBtwH.AddDat(TIntPr(NI.GetId(), NI.GetOutNId(e)), 0);
        } else {
          // add each edge only once in undirected graphs
          if (NI.GetId() < NI.GetOutNId(e)) {
            EdgeBtwH.AddDat(TIntPr(NI.GetId(), NI.GetOutNId(e)), 0); 
          }
        }
      }
      // add incoming edges in directed graphs that were not added yet
      if (Graph->HasFlag(gfDirected) && !IsDir) {
        for (int e = 0; e < NI.GetInDeg(); e++) {
          if (NI.GetId() < NI.GetInNId(e)  &&
              !Graph->IsEdge(NI.GetId(), NI.GetInNId(e))) {
            EdgeBtwH.AddDat(TIntPr(NI.GetId(), NI.GetInNId(e)), 0);  
          } 
        }
      }
    }
    sigma.AddDat(NI.GetId(), 0);
    d.AddDat(NI.GetId(), -1);
    P.AddDat(NI.GetId(), TIntV());
    delta.AddDat(NI.GetId(), 0);
  }
  // calc betweeness
  for (int k=0; k < BtwNIdV.Len(); k++) {
    const PNEANet::TObj::TNodeI NI = Graph->GetNI(BtwNIdV[k]);
    // reset
    for (int i = 0; i < sigma.Len(); i++) {
      sigma[i]=0;  d[i]=-1;  delta[i]=0;  P[i].Clr(false);
    }
    S.Clr(false);
    Q.Clr(false);
    sigma.AddDat(NI.GetId(), 1);
    d.AddDat(NI.GetId(), 0);
    Q.Push(NI.GetId());
    while (! Q.Empty()) {
      const int v = Q.Top();  Q.Pop();
      const PNEANet::TObj::TNodeI NI2 = Graph->GetNI(v);
      S.Push(v);
      const double VDat = d.GetDat(v);
      // iterate over all outgoing edges
      for (int e = 0; e < NI2.GetOutDeg(); e++) {
        const int w = NI2.GetOutNId(e);
        const int eid = NI2.GetOutEId(e);

        if (d.GetDat(w) < 0) { // find w for the first time
          Q.Push(w);
          d.AddDat(w, VDat+Attr[eid]);
        }
        //shortest path to w via v ?
        if (d.GetDat(w) == VDat+Attr[eid]) {
          sigma.AddDat(w) += sigma.GetDat(v);
          P.GetDat(w).Add(v);
        }
      }
      // if ignoring direction in directed networks, iterate over incoming edges
      if (Graph->HasFlag(gfDirected) && !IsDir) {
        for (int e = 0; e < NI2.GetInDeg(); e++) {
          const int w = NI2.GetInNId(e);
          // skip neighbors that are also outgoing
          if (Graph->IsEdge(NI2.GetId(), w)) {
            continue;
          }
          const int eid = NI2.GetInEId(e);

          if (d.GetDat(w) < 0) { // find w for the first time
            Q.Push(w);
            d.AddDat(w, VDat+Attr[eid]);
          }
          //shortest path to w via v ?
          if (d.GetDat(w) == VDat+Attr[eid]) {
            sigma.AddDat(w) += sigma.GetDat(v);
            P.GetDat(w).Add(v);
          }
        }
      }
    }
    
    while (! S.Empty()) {
      const int w = S.Top();
      const double SigmaW = sigma.GetDat(w);
      const double DeltaW = delta.GetDat(w);
      const TIntV NIdV = P.GetDat(w);
      S.Pop();
      for (int i = 0; i < NIdV.Len(); i++) {
        const int NId = NIdV[i];
        const double c = (sigma.GetDat(NId)*1.0/SigmaW) * (1+DeltaW);
        delta.AddDat(NId) += c;
        if (DoEdgeCent) {
          if (Graph->HasFlag(gfDirected) && IsDir) {
            EdgeBtwH.AddDat(TIntPr(NId, w)) += c;
          } else {
            EdgeBtwH.AddDat(TIntPr(TMath::Mn(NId, w), TMath::Mx(NId, w))) += c;
          }
        }
      }
      if (DoNodeCent && w != NI.GetId()) {
        NodeBtwH.AddDat(w) += delta.GetDat(w)/2.0; }
    }
  }
}

void GetWeightedBetweennessCentr(const PNEANet Graph, TIntFltH& NodeBtwH, TIntPrFltH& EdgeBtwH, const TFltV& Attr, const double& NodeFrac, const bool& IsDir) {
  TIntV NIdV;  Graph->GetNIdV(NIdV);
  if (NodeFrac < 1.0) { // calculate beetweenness centrality for a subset of nodes
    NIdV.Shuffle(TInt::Rnd);
    for (int i = int((1.0-NodeFrac)*NIdV.Len()); i > 0; i--) {
      NIdV.DelLast(); }
  }
  GetWeightedBetweennessCentr(Graph, NIdV, NodeBtwH, true, EdgeBtwH, true,
    Attr, IsDir);
}

void GetWeightedBetweennessCentr(const PNEANet Graph, TIntFltH& NodeBtwH, const TFltV& Attr, const double& NodeFrac, const bool& IsDir) {
  TIntPrFltH EdgeBtwH;
  TIntV NIdV;  Graph->GetNIdV(NIdV);
  if (NodeFrac < 1.0) { // calculate beetweenness centrality for a subset of nodes
    NIdV.Shuffle(TInt::Rnd);
    for (int i = int((1.0-NodeFrac)*NIdV.Len()); i > 0; i--) {
      NIdV.DelLast(); }
  }
  GetWeightedBetweennessCentr(Graph, NIdV, NodeBtwH, true, EdgeBtwH, false,
    Attr, IsDir);
}

void GetWeightedBetweennessCentr(const PNEANet Graph, TIntPrFltH& EdgeBtwH, const TFltV& Attr, const double& NodeFrac, const bool& IsDir) {
  TIntFltH NodeBtwH;
  TIntV NIdV;  Graph->GetNIdV(NIdV);
  if (NodeFrac < 1.0) { // calculate beetweenness centrality for a subset of nodes
    NIdV.Shuffle(TInt::Rnd);
    for (int i = int((1.0-NodeFrac)*NIdV.Len()); i > 0; i--) {
      NIdV.DelLast(); }
  }
  GetWeightedBetweennessCentr(Graph, NIdV, NodeBtwH, false, EdgeBtwH, true,
    Attr, IsDir);
}

/// Gets sequence of PageRank tables from given \c GraphSeq.
TTableIterator GetMapPageRank(
    const TVec<PNEANet>& GraphSeq,
    TTableContext* Context,
    const double& C = 0.85, const double& Eps = 1e-4, const int& MaxIter = 100) {
  TVec<PTable> TableSeq(GraphSeq.Len());
  TSnap::MapPageRank(GraphSeq, TableSeq, Context, C, Eps, MaxIter);
  return TTableIterator(TableSeq);
}

/// Gets sequence of Hits tables from given \c GraphSeq.
TTableIterator GetMapHitsIterator(
    const TVec<PNEANet>& GraphSeq,
    TTableContext* Context,
    const int& MaxIter = 20) {
  TVec<PTable> TableSeq(GraphSeq.Len());
  TSnap::MapHits(GraphSeq, TableSeq, Context, MaxIter);
  return TTableIterator(TableSeq);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//added by wangjufan
  
struct TFrontierMN {
    double SrcDist;
    unsigned long SrcNID;
//    unsigned long PrtNID;
};
  
int GetWeightedShortestPathByDijkstraMemoryArrayHash (IShortestPathGraph* Graph,
                                                      unsigned long startNId,
                                             TIntFltH& NIdDistH,
                                             TDijkstraStat& stat,
                                                std::unordered_map<unsigned long, unsigned long>& NodeID2CNodeID) {
    std::vector<TFrontierMN> frontier;
    auto comparator = [&] (TFrontierMN& left, TFrontierMN& right) {
        return left.SrcDist > right.SrcDist;
    };
  
  BitsContainerType *NodeDistFlag = Graph->getDistFlagVector();//opt
  
//    unsigned long maxnodeid = Graph->getMaxNodeID()+2;
//   double * NodeDistH = ( double *)malloc(sizeof( double)*maxnodeid);
//  for (unsigned long i=0; i < maxnodeid; i++) {
//    NodeDistH[i] = __DBL_MAX__;
//  }

    std::unordered_map<int, double> NodeDistH = Graph->getNodeDistH();
  
    NodeDistH[startNId] = 0;
    struct TFrontierMN node = {0, startNId};
    frontier.push_back(node);
    
    while (! frontier.empty()) {
      
      if (frontier.size() > stat.getUseless()) {
        stat.incrUseless(frontier.size());
      }
      
        pop_heap(begin(frontier), end(frontier), comparator);
        auto& frontierNode = frontier.back();
        frontier.pop_back();
        int&& SrcNID = frontierNode.SrcNID;
        double ParentDistance = frontierNode.SrcDist;
        double& plen = NodeDistH[SrcNID];//45
        if (plen != ParentDistance) {
            continue;
        }
      
//      unsigned long ClusteringNID = frontierNode.SrcNID;///opt
//      unsigned long pos = ClusteringNID & 0x3F;
//      unsigned long index = ClusteringNID >> 6;
//      SET_BIT(NodeDistFlag[index], pos);
                
        std::vector<TEdgeTuple>& AttrFltIntKV =  Graph->GetSortedAttrByNode(SrcNID);
 
//      std::for_each(AttrFltIntKV.begin(), AttrFltIntKV.end(), [&](TEdgeTuple& tpl) {
//        long int && dd =tpl.to;
//        double& plen = NodeDistH[dd];
//        double DstDistance = ParentDistance + tpl.EdgeLen;//origin
//        stat.incrEdgeNum();
//          if (plen > DstDistance) [[unlikely]] {
//              NodeDistH[tpl.to] = DstDistance;
//              struct TFrontierMN node = {DstDistance, tpl.to};
//              frontier.push_back(node);
//              push_heap(begin(frontier), end(frontier), comparator);
//            stat.incrHeapItem();//3
//          }
//          });
      
      ///////
      int currentEdge = 0;
      unsigned long size = AttrFltIntKV.size();
      while (size >  currentEdge) {
        TEdgeTuple& tpl = (AttrFltIntKV)[currentEdge];
        unsigned long&     dd =tpl.to;
        double&  plen = NodeDistH[dd];
        double DstDistance = ParentDistance + tpl.EdgeLen;//origin
        currentEdge+=1;
        stat.incrEdgeNum();
          if (plen > DstDistance) [[unlikely]]
          {
              NodeDistH[dd] = DstDistance;
//              struct TFrontierMN node = {DstDistance, dd};
//              frontier.push_back(node);
//              push_heap(begin(frontier), end(frontier), comparator);
//            stat.incrHeapItem();//3
          }
        struct TFrontierMN node = {DstDistance, dd};
        frontier.push_back(node);
        push_heap(begin(frontier), end(frontier), comparator);
      stat.incrHeapItem();//3
      }
      
      /////////
//      if(size >  currentEdge) {
//        TEdgeTuple tpl = (AttrFltIntKV)[currentEdge];
//        unsigned long     dd =tpl.to;
//        double  plen = NodeDistH[dd];
//        double DstDistance = ParentDistance + tpl.EdgeLen;//origin
//
//        TEdgeTuple tplNext ;
//        unsigned long  ddNext ;
//        double plenNext ;
//        double DstDistanceNext ;//origin
//        do {
//          currentEdge+=1;
//          if (size >  currentEdge) {
//            tplNext = (AttrFltIntKV)[currentEdge];
//            ddNext =tplNext.to;
//            plenNext = NodeDistH[ddNext];
//            DstDistanceNext = ParentDistance + tplNext.EdgeLen;//origin
//          }
//          stat.incrEdgeNum();
//            if (plen > DstDistance) [[unlikely]] {
//                NodeDistH[dd] = DstDistance;
//                struct TFrontierMN node = {DstDistance, dd};
//                frontier.push_back(node);
//                push_heap(begin(frontier), end(frontier), comparator);
//              stat.incrHeapItem();//3
//            }
//             dd = ddNext;
//           plen = plenNext;
//           DstDistance = DstDistanceNext;
//        }while (size >  currentEdge);
//      }
      
  }
    
    NIdDistH.Clr(false);
//  for (unsigned long i = 0; i < maxnodeid; i++) {
//      if (NodeDistH[i] != __DBL_MAX__) {
//          NIdDistH.AddDat(i, NodeDistH[i]);
//      }
//  }
    for (auto x : NodeDistH) {
        if (x.second != __DBL_MAX__) {
//            stat.incrNodeCount();
            NIdDistH.AddDat(x.first, x.second);
        }
    }//wjf
    return 0;
}


/////////////////////
struct TFrontierSSNode {
    double DstDist;
    unsigned long DstNID;
    unsigned long PrtNID;
    double PrtDist;
    unsigned long cto;
    int EdgeIndex;
    std::vector<TEdgeTuple>* AttrFltIntKVPtr;
    int size;
};

struct TDistanceNode {
    unsigned long NID;
    double Dist;
};

//inline void pushNextFrontierNodeBest(std::vector<TFrontierSSNode>&  frontier,
//                                     TSmallStepStat& stat,
//                                     struct TFrontierSSNode& fnode,
//                                     BitsContainerType* NodeDistFlag) {
//    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
//        return left.DstDist > right.DstDist;
//    };
//    int const Size = fnode.size;
//    while (fnode.EdgeIndex < Size) {
//        TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[fnode.EdgeIndex++];
//        BitsContainerType vv = NodeDistFlag[tpl.index()];
//      fnode.DstNID = tpl.to;
//      bool flag = GET_BIT(vv, tpl.pos());
//        if (!flag)  [[unlikely]]{
//          fnode.DstDist = tpl.EdgeLen + fnode.PrtDist;
//            frontier.push_back(fnode);
//            push_heap(begin(frontier), end(frontier), comparator);
//            break;
//        }
//    }
//}
//
//int GetWeightedShortestPathBySmallStepOnNGraphBest(IShortestPathGraph* Graph,
//                                                   unsigned long SrcNId,
//                                               TIntFltH& NIdDistH,
//                                               TSmallStepStat& stat) {
//    NIdDistH.Clr(false);
//    std::vector<TFrontierSSNode> frontier;
//    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
//        return left.DstDist > right.DstDist;
//    };
//    BitsContainerType *NodeDistFlag = Graph->getDistFlagVector();
//
//    unsigned long index = SrcNId >> BitsShiftCount;
//    unsigned long pos = SrcNId & BitsMaskForEquivalenceClasses;
//    SET_BIT(NodeDistFlag[index], pos);
//    NIdDistH.AddDat(SrcNId, 0);
//
//    std::vector<TDistanceNode> distanceVector;
//
//    struct TFrontierSSNode pnode;
//    pnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(SrcNId);
//    pnode.size = pnode.AttrFltIntKVPtr->size();
//    pnode.EdgeIndex = 0;
//    pnode.PrtDist = 0;
//    pnode.PrtNID = 0;
//    TSnap::pushNextFrontierNodeBest(frontier, stat, pnode, NodeDistFlag);
//
//    while (!frontier.empty()) {
//        pop_heap(begin(frontier), end(frontier), comparator);
//        auto frontierNode = frontier.back();
//        frontier.pop_back();
//        unsigned long DstNID = frontierNode.DstNID;
//        unsigned long pos = DstNID & BitsMaskForEquivalenceClasses;
//        unsigned long index = DstNID >> BitsShiftCount;
//
//        if (! GET_BIT(NodeDistFlag[index], pos))  [[unlikely]]{
//            SET_BIT(NodeDistFlag[index], pos);
//            stat.incrNodeCount();
//            distanceVector.push_back({DstNID,frontierNode.DstDist});
//            struct TFrontierSSNode nnode;
//            nnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(frontierNode.DstNID);
//            nnode.size = nnode.AttrFltIntKVPtr->size();
//            nnode.EdgeIndex = 0;
//            nnode.PrtDist = frontierNode.DstDist;
//            nnode.PrtNID = frontierNode.DstNID;
//            TSnap::pushNextFrontierNodeBest(frontier, stat, nnode, NodeDistFlag);
//        }
////      stat.incrInvalidHeapItem();//1
//        TSnap::pushNextFrontierNodeBest(frontier, stat, frontierNode, NodeDistFlag);
//    }
//    for (auto i = distanceVector.cbegin(); i != distanceVector.cend(); ++i) {
//        NIdDistH.AddDat(i->NID, i->Dist);
//    }//wjf
//    return 0;
//}

static struct TFrontierSSNode preNode;
////////////////////////////////////////////
inline void pushNextFrontierNodeBestWithClustering(std::vector<TFrontierSSNode>&  frontier,
                                     TSmallStepStat& stat,
                                     struct TFrontierSSNode& fnode,
                                     BitsContainerType* NodeDistFlag,
                        std::unordered_map<unsigned long, unsigned long>& NodeID2CNodeID) {
    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
        return left.DstDist > right.DstDist;
    };
    int const Size = fnode.size;
  
  unsigned long ClusteringNID;
  unsigned long index;
  unsigned long pos;
  
  if (fnode.EdgeIndex < Size) {
    TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[fnode.EdgeIndex];
     ClusteringNID = tpl.cto;
     index = ClusteringNID >> BitsShiftCount;
     pos = ClusteringNID & BitsMaskForEquivalenceClasses;
    while (fnode.EdgeIndex < Size) {
      TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[fnode.EdgeIndex];
      fnode.EdgeIndex++;
      BitsContainerType& vv = NodeDistFlag[index];
      if (!GET_BIT(vv, pos)) [[unlikely]]
      {
        fnode.DstNID = tpl.to;
        fnode.cto = tpl.cto;
        fnode.DstDist = tpl.EdgeLen + fnode.PrtDist;
        frontier.push_back(fnode);
        stat.incrHeapItem();//2
        push_heap(begin(frontier), end(frontier), comparator);
        break;
      } else {
        fnode.DstNID = tpl.to;
        fnode.cto = tpl.cto;
        fnode.DstDist = tpl.EdgeLen + fnode.PrtDist;
        frontier.push_back(fnode);
        stat.incrHeapItem();//2
        push_heap(begin(frontier), end(frontier), comparator);
        break;
//        TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[fnode.EdgeIndex];
//         ClusteringNID = tpl.cto;
//         index = ClusteringNID >> BitsShiftCount;
//         pos = ClusteringNID & BitsMaskForEquivalenceClasses;
      }
    }
  }
}

int GetWeightedShortestPathBySmallStepOnNGraphBestWithClustering(IShortestPathGraph* Graph,
                                                                 unsigned long SrcNId,
                                               TIntFltH& NIdDistH,
                                               TSmallStepStat& stat,
                                  std::unordered_map<unsigned long, unsigned long>& NodeID2CNodeID) {
    NIdDistH.Clr(false);
    std::vector<TFrontierSSNode> frontier;
    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
        return left.DstDist > right.DstDist;
    };
    BitsContainerType *NodeDistFlag = Graph->getDistFlagVector();
    
    unsigned long ClusteringNID = NodeID2CNodeID[SrcNId];
    unsigned long index = ClusteringNID >> BitsShiftCount;
    unsigned long pos = ClusteringNID & BitsMaskForEquivalenceClasses;
    SET_BIT(NodeDistFlag[index], pos);
    NIdDistH.AddDat(SrcNId, 0);
    stat.incrNodeCount();
  
    std::vector<TDistanceNode> distanceVector;
    
    struct TFrontierSSNode pnode;
    pnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(SrcNId);
    pnode.size = pnode.AttrFltIntKVPtr->size();
    pnode.EdgeIndex = 0;
    pnode.PrtDist = 0;
    pnode.PrtNID = 0;
    TSnap::pushNextFrontierNodeBestWithClustering(frontier, stat, pnode, NodeDistFlag, NodeID2CNodeID);
   
  double PathLen = 0;
    while (!frontier.empty()) {
      
      if (frontier.size() > stat.getUseless()) {
        stat.incrUseless(frontier.size());
      }
      
        pop_heap(begin(frontier), end(frontier), comparator);
        auto frontierNode = frontier.back();
        frontier.pop_back();
        unsigned long DstNID = frontierNode.DstNID;
      
      unsigned long ClusteringNID = frontierNode.cto;
      assert(ClusteringNID == frontierNode.cto);
        unsigned long pos = ClusteringNID & BitsMaskForEquivalenceClasses;
        unsigned long index = ClusteringNID >> BitsShiftCount;
      
        if (! GET_BIT(NodeDistFlag[index], pos)) [[unlikely]]{
          PathLen = frontierNode.DstDist;
            SET_BIT(NodeDistFlag[index], pos);
            stat.incrNodeCount();
            distanceVector.push_back({DstNID,frontierNode.DstDist});
            struct TFrontierSSNode nnode;
            nnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(frontierNode.DstNID);
            nnode.size = nnode.AttrFltIntKVPtr->size();
            nnode.EdgeIndex = 0;
            nnode.PrtDist = frontierNode.DstDist;
            nnode.PrtNID = frontierNode.DstNID;
            TSnap::pushNextFrontierNodeBestWithClustering(frontier, stat, nnode, NodeDistFlag, NodeID2CNodeID);
          
          TSnap::pushNextFrontierNodeBestWithClustering(frontier, stat, frontierNode, NodeDistFlag, NodeID2CNodeID);
          
        }else {
          TSnap::pushNextFrontierNodeBestWithClustering(frontier, stat, frontierNode, NodeDistFlag, NodeID2CNodeID);
        }
        
    }
    for (auto i = distanceVector.cbegin(); i != distanceVector.cend(); ++i) {
        NIdDistH.AddDat(i->NID, i->Dist);
    }//wjf
    return 0;
}


/////////////////////////////////////////// id to c id
inline void pushNextFrontierNode4Clustering(std::vector<TFrontierSSNode>&  frontier,
                                     struct TFrontierSSNode& fnode,
                                     BitsContainerType* NodeDistFlag,
                                            std::vector<unsigned long>& WaitingNodeIDs) {
    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
        return left.DstDist > right.DstDist;
    };
    int const Size = fnode.size;
  
  if(fnode.EdgeIndex == 0) {
    int i = 0;
    while (i < Size) {
        TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[i++];
        BitsContainerType vv = NodeDistFlag[tpl.index()];
        bool flag = GET_BIT(vv, tpl.pos());
        if (!flag) {
          WaitingNodeIDs.push_back(tpl.to);
        }
    }
  }
  
    while (fnode.EdgeIndex < Size) {
        TEdgeTuple& tpl = (*fnode.AttrFltIntKVPtr)[fnode.EdgeIndex++];
        BitsContainerType vv = NodeDistFlag[tpl.index()];
      fnode.DstNID = tpl.to;
      bool flag = GET_BIT(vv, tpl.pos());
        if (!flag)  [[unlikely]]{
            fnode.DstDist = tpl.EdgeLen + fnode.PrtDist;
            frontier.push_back(fnode);
            push_heap(begin(frontier), end(frontier), comparator);
            break;
        }
    }
}
std::vector<unsigned long> GetWeightedShortestPathPoints(
                                  IShortestPathGraph* Graph,
                                  unsigned long SrcNId,
                                  BitsContainerType *NodeDistFlag,
                                  std::vector<unsigned long>& WaitingNodeIDs) {
  std::vector<unsigned long> members;//待选等价类
  
    std::vector<TFrontierSSNode> frontier;
    auto comparator = [&] (TFrontierSSNode& left, TFrontierSSNode& right) {
        return left.DstDist > right.DstDist;
    };
    
    unsigned long index = SrcNId >> BitsShiftCount;
    unsigned long pos = SrcNId & BitsMaskForEquivalenceClasses;
    SET_BIT(NodeDistFlag[index], pos);
    members.push_back(SrcNId);
  
    std::vector<TDistanceNode> distanceVector;//最短路径长度
    
    struct TFrontierSSNode pnode;
    pnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(SrcNId);
    pnode.size = pnode.AttrFltIntKVPtr->size();
    pnode.EdgeIndex = 0;
    pnode.PrtDist = 0;
    pnode.PrtNID = 0;
    TSnap::pushNextFrontierNode4Clustering(frontier, pnode, NodeDistFlag, WaitingNodeIDs);
   
    while (!frontier.empty() & members.size() < EquivalenceClassesElementNum) {
        pop_heap(begin(frontier), end(frontier), comparator);
        auto frontierNode = frontier.back();
        frontier.pop_back();
        unsigned long DstNID = frontierNode.DstNID;
        unsigned long pos = DstNID & BitsMaskForEquivalenceClasses;
        unsigned long index = DstNID >> BitsShiftCount;
      
        if (! GET_BIT(NodeDistFlag[index], pos))  [[unlikely]]{
          SET_BIT(NodeDistFlag[index], pos);
          members.push_back(DstNID);
            distanceVector.push_back({DstNID,frontierNode.DstDist});
            struct TFrontierSSNode nnode;
            nnode.AttrFltIntKVPtr =  &Graph->GetSortedAttrByNode(frontierNode.DstNID);
            nnode.size = nnode.AttrFltIntKVPtr->size();
            nnode.EdgeIndex = 0;
            nnode.PrtDist = frontierNode.DstDist;
            nnode.PrtNID = frontierNode.DstNID;
            TSnap::pushNextFrontierNode4Clustering(frontier, nnode, NodeDistFlag, WaitingNodeIDs);
        }
        TSnap::pushNextFrontierNode4Clustering(frontier, frontierNode, NodeDistFlag, WaitingNodeIDs);
    }
    return members;
}

void nodeID2ClusteringID(std::unordered_map<unsigned long, unsigned long>& nid2clustering,
                         unsigned long SrcNId,
                         IShortestPathGraph* Graph,
                         BitsContainerType *NodeDistFlag,
                         std::vector<std::vector<unsigned long>>& equalClss) {
  //待选节点 K=2 每次最多选择两个节点
  //MarkNum = 64 计算最短路径64个节点的最短路径
  std::vector<unsigned long> WaitingNodeIDs;   //待选 迭代起点
  WaitingNodeIDs.push_back(SrcNId);
  
  while (WaitingNodeIDs.size()) {
    unsigned long rootID = WaitingNodeIDs.back();
    WaitingNodeIDs.pop_back();
    unsigned long pos = rootID & BitsMaskForEquivalenceClasses;
    unsigned long index = rootID >> BitsShiftCount;
    if (GET_BIT(NodeDistFlag[index], pos)) {
      continue;
    }
    auto EC = GetWeightedShortestPathPoints(Graph, rootID, NodeDistFlag, WaitingNodeIDs);
    equalClss.push_back(EC);
  }
}

/////////////////////
//动态生产图

}; // namespace TSnap
 
