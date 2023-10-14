// Small example testing basic functionality of SNAP

#include "Snap.h"

#include <time.h>
#include <cstdlib>
#include <fstream>

#include <streambuf>
#include <cmath>

#include <thread>
#include <stdlib.h>

void testHeap( IShortestPathGraph * Graph1, int nid,
          clock_t* time1 ,clock_t* time2,
          TSnap::TDijkstraStat& DijkstraStat,
          TSnap::TSmallStepStat& SmallStepStat) {
//    Graph1->getMaxNodeID();
//    TIntFltH NIdDistDj;
//    clock_t startTime2,endTime2;
//    startTime2 = clock();
//    TSnap::GetWeightedShortestPathByDijkstraMemoryArrayHash(Graph1, nid, NIdDistDj,
//                                                   DijkstraStat);
//
//
//    endTime2 = clock();
//    *time2 = (double)(endTime2 - startTime2) ;
//    TIntFltH NIdDistStep;
//   clock_t startTime,endTime;
//       startTime = clock();
//   TSnap::GetWeightedShortestPathBySmallStepOnNGraphBest(Graph1, nid, NIdDistStep,
//                                                     SmallStepStat);
//       endTime = clock();
//    *time1 = (double)(endTime - startTime);
//    for (TIntFltH::TIter It = NIdDistStep.BegI(); It < NIdDistStep.EndI(); It++) {
//      int node_id = It.GetKey();//node id
//      double centrStep = It.GetDat();//distanc to source node
//        if (centrStep != __DBL_MAX__) {
//            double centrDj = NIdDistDj.GetDat(node_id);
//            assert(centrStep == centrDj);
//        }
//    }
//    for (TIntFltH::TIter It = NIdDistDj.BegI(); It < NIdDistDj.EndI(); It++) {
//      int node_id = It.GetKey();//node id
//      double centrDj = It.GetDat();//distanc to source node
//        if (centrDj != __DBL_MAX__) {
//            double centrStep = NIdDistStep.GetDat(node_id);
//            assert(centrDj == centrStep);
//        }
//    }
}

void compairSmallStepAndDijkstras(std::string name, bool dir) {

//    std::string bstr = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/u-w/";
//    std::string cstr = name + "/" + name;
//
//    std::string dbstr = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/d-w/";
//    std::string dcstr = name;
//
//    std::string istr = ".edges";
//
//
//    PNGraph Graph1;
//    PUNGraph UGraph1;
//    if (dir) {
//        Graph1 = TSnap::LoadAttrEdgeList<PNGraph>((dbstr+dcstr+istr).c_str(), 0, 1, 2);
//        Graph1->sortEdgeByAttr();
//    }else {
//        UGraph1 = TSnap::LoadAttrEdgeList<PUNGraph>((bstr+cstr+istr).c_str(), 0, 1, 2);
//        UGraph1->sortEdgeByAttr();
//    }
//
//
//        srand((unsigned)time(NULL));
//      int *randomNodeIDIndex;
//    if (dir) {
//        randomNodeIDIndex= (int*)malloc(sizeof(int)*Graph1->GetNodes());
//        memset(randomNodeIDIndex,0, sizeof(int)*Graph1->GetNodes());
//    }else {
//        randomNodeIDIndex= (int*)malloc(sizeof(int)*UGraph1->GetNodes());
//        memset(randomNodeIDIndex,0, sizeof(int)*UGraph1->GetNodes());
//    }
//
//        int count = 0;
//       int account = 0;
//    int rindex ;
//        while (count < 10 && account < 2000) {
//            if (dir) {
//                rindex = (rand()%Graph1->GetNodes());
//            }else {
//                rindex = (rand()%UGraph1->GetNodes());
//            }
//            if (randomNodeIDIndex[rindex] ==0) {
//                randomNodeIDIndex[rindex] = 1;
//                count++;
//            }
//            account++;
//        }
//
//        clock_t smallStepTime = 0;
//        clock_t DijkstraTime = 0;
//        TSnap::TDijkstraStat DijkstraStat;
//        TSnap::TSmallStepStat SmallStepStat;
//        std::ofstream oFile;
//    std::string genbase = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/gen/";
//    std::string ostr = ".csv";
//    oFile.open((genbase+name+ostr).c_str(), std::ios::out | std::ios::trunc);
//
//  oFile<< "节点" << ","
//       << "S访问边" << ","
//       << "D ReHeap" << ","
//       << "nid"<< ","
//       <<  "D/S"<< ","
//        << "S-CLOCKS" << "," << "S-time" << ","
//        << "D-CLOCKS" << "," << "D-time" << ","
//        << "实现方式"
//       << std::endl;
//
//    int index = 0;
//    if (dir) {
//        for (TNGraph::TNodeI NI = Graph1->BegNI(); NI < Graph1->EndNI(); NI++) {
//            if(randomNodeIDIndex[index] == 1){
//                int nid = NI.GetId();
//                DijkstraStat.restat();
//                SmallStepStat.restat();
//
//                DijkstraStat.restat();
//                SmallStepStat.restat();
//                auto graph = Graph1();
//                testHeap(graph, nid,
//                     &smallStepTime, &DijkstraTime,
//                     DijkstraStat, SmallStepStat);
//              oFile<< DijkstraStat.getNodeCount() << ","
//                << (double)SmallStepStat.getVisitedEdgeNum() << ","
//                << (double)DijkstraStat.getReHeapCount() << ","
//                << nid << ","
//                << (double)((double)DijkstraTime/(double)smallStepTime)<< ","
//                << smallStepTime << "," << (double)smallStepTime/CLOCKS_PER_SEC << ","
//                << DijkstraTime << "," << (double)DijkstraTime/CLOCKS_PER_SEC <<","
//                << "heap"
//                << std::endl;
//
//            }
//            index++;
//        }
//    }else {
//        for (TUNGraph::TNodeI NI = UGraph1->BegNI(); NI < UGraph1->EndNI(); NI++) {
//            if(randomNodeIDIndex[index] == 1){
//                int nid = NI.GetId();
//                DijkstraStat.restat();
//                SmallStepStat.restat();
//
//                DijkstraStat.restat();
//                SmallStepStat.restat();
//                testHeap(UGraph1(), nid,
//                     &smallStepTime, &DijkstraTime,
//                     DijkstraStat, SmallStepStat);
//           oFile << DijkstraStat.getNodeCount() << ","
//                << (double)SmallStepStat.getVisitedEdgeNum()/DijkstraStat.getNodeCount() << ","
//                << (double)DijkstraStat.getReHeapCount()/DijkstraStat.getNodeCount() << ","
//                << nid << ","
//              << (double)((double)DijkstraTime/(double)smallStepTime)<< ","
//                << smallStepTime << "," << (double)smallStepTime/CLOCKS_PER_SEC << ","
//                << DijkstraTime << "," << (double)DijkstraTime/CLOCKS_PER_SEC <<","
//                << "heap"
//                << std::endl;
//
//            }
//            index++;
//
//        }
//    }
//
//    oFile.close();
}
//////////////////////////////////////////////
void testHeapWithClustering(IShortestPathGraph * Graph1, int nid,
          clock_t* time1 ,clock_t* time2,
          TSnap::TDijkstraStat& DijkstraStat,
          TSnap::TSmallStepStat& SmallStepStat,
          std::unordered_map<unsigned long, unsigned long>& NodeID2CNodeID) {
  
    Graph1->getMaxNodeID();
    TIntFltH NIdDistDj;
    clock_t startTime2,endTime2;
    startTime2 = clock();
    TSnap::GetWeightedShortestPathByDijkstraMemoryArrayHash(Graph1, nid, NIdDistDj,
                                                   DijkstraStat,NodeID2CNodeID);
    
    
    endTime2 = clock();
    *time2 = (double)(endTime2 - startTime2) ;
    TIntFltH NIdDistStep;
   clock_t startTime,endTime;
       startTime = clock();
   TSnap::GetWeightedShortestPathBySmallStepOnNGraphBestWithClustering(Graph1, nid, NIdDistStep,
                                                     SmallStepStat, NodeID2CNodeID);
       endTime = clock();
    *time1 = (double)(endTime - startTime);
  int diffcount = 0;
    for (TIntFltH::TIter It = NIdDistStep.BegI(); It < NIdDistStep.EndI(); It++) {
      int node_id = It.GetKey();//node id
      double centrStep = It.GetDat();//distanc to source node
      double centrDj = NIdDistDj.GetDatWithDefault(node_id, -1);
      if (centrStep != centrDj) {
        diffcount++;
      }
//      assert(centrStep == centrDj);
    }
    for (TIntFltH::TIter It = NIdDistDj.BegI(); It < NIdDistDj.EndI(); It++) {
      int node_id = It.GetKey();//node id
      double centrDj = It.GetDat();//distanc to source node
      double centrStep = NIdDistStep.GetDatWithDefault(node_id, -1);
      
      //            assert(centrDj == centrStep);
                  if (centrStep != centrDj) {
                    diffcount++;
                  }
    }
  if (diffcount > 0) {
    printf("\n===================\n", diffcount);
    printf("diff count = %d", diffcount);
    printf("\n===================\n", diffcount);
  }
}
void produceClusteringID(std::string name, bool dir) {
//  std::string bstr = "/Volumes/build/";
  std::string bstr = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/u-w/";
  std::string cstr = name + "/" + name;
  
//  std::string dbstr = "/Volumes/build/";
  std::string dbstr = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/d-w/";
  std::string dcstr = name;
  
  std::string istr = ".edges";
  
  
  PNGraph Graph1;
  PUNGraph UGraph1;
  if (dir) {
      Graph1 = TSnap::LoadAttrEdgeList<PNGraph>((dbstr+dcstr+istr).c_str(), 0, 1, 2);
      Graph1->sortEdgeByAttr();
  }else {
      UGraph1 = TSnap::LoadAttrEdgeList<PUNGraph>((bstr+cstr+istr).c_str(), 0, 1, 2);
      UGraph1->sortEdgeByAttr();
  }
      

      srand((unsigned)time(NULL));
    int *randomNodeIDIndex;
  int nodeNUm = 0;
  if (dir) {
    nodeNUm = Graph1->GetNodes();
      randomNodeIDIndex= (int*)malloc(sizeof(int)*Graph1->GetNodes());
      memset(randomNodeIDIndex,0, sizeof(int)*Graph1->GetNodes());
  }else {
    nodeNUm = UGraph1->GetNodes();
      randomNodeIDIndex= (int*)malloc(sizeof(int)*UGraph1->GetNodes());
      memset(randomNodeIDIndex,0, sizeof(int)*UGraph1->GetNodes());
  }
      
      int count = 0;
     int account = 0;
  int rindex ;
      while (count < 10 && account < 2000) {
          if (dir) {
              rindex = (rand()%Graph1->GetNodes());
          }else {
              rindex = (rand()%UGraph1->GetNodes());
          }
          if (randomNodeIDIndex[rindex] ==0) {
              randomNodeIDIndex[rindex] = 1;
              count++;
          }
          account++;
      }

  clock_t smallStepTime = 0;
  clock_t DijkstraTime = 0;
  TSnap::TDijkstraStat DijkstraStat;
  TSnap::TSmallStepStat SmallStepStat;

  std::ofstream clusteringFile;
std::string genbase = "/Users/jufanwang/SmallStepShortestPath/snap/dataset/gen/";
  std::string ostr = ".csv";
  std::string clstr = "-clustering";
  clusteringFile.open((genbase+name+clstr+ostr).c_str(), std::ios::out | std::ios::trunc);

  clusteringFile<< "节点" << ","<< "C节点" << std::endl;
  
  std::unordered_map<unsigned long, unsigned long> NodeID2CNodeID;
  
  BitsContainerType *NodeDistFlag = NULL;
  if (dir) {
    NodeDistFlag = Graph1->getDistFlagVector();//标记
  }else {
    NodeDistFlag = UGraph1->getDistFlagVector();//标记
  }
  
  int index = 0;
  std::vector<std::vector<unsigned long>> equalClss;  //等价类
  if (dir) {
    for (TNGraph::TNodeI NI = Graph1->BegNI(); NI < Graph1->EndNI(); NI++) {
      int nid = NI.GetId();
      unsigned long pos = nid & BitsMaskForEquivalenceClasses;
      unsigned long index = nid >> BitsShiftCount;
      if (GET_BIT(NodeDistFlag[index], pos)) {
        continue;
      }
      TSnap::nodeID2ClusteringID(NodeID2CNodeID, nid, Graph1(), NodeDistFlag, equalClss);
    }
  }else {
    for (TUNGraph::TNodeI NI = UGraph1->BegNI(); NI < UGraph1->EndNI(); NI++) {
      int nid = NI.GetId();
      unsigned long pos = nid & BitsMaskForEquivalenceClasses;
      unsigned long index = nid >> BitsShiftCount;
      if (GET_BIT(NodeDistFlag[index], pos)) {
        continue;
      }
      TSnap::nodeID2ClusteringID(NodeID2CNodeID, nid, UGraph1(), NodeDistFlag, equalClss);
    }
  }
  
  //编码 i = 0 等价类 j = 0-63类内元素，转成ID
//  std::vector<unsigned long> sum;
  std::vector<std::vector<unsigned long>> smallEqualClss;
  unsigned long clssID = 0;
  for (std::vector<unsigned long> clss : equalClss) {
    int size = clss.size();
    if (size == EquivalenceClassesElementNum) {
      int eindex = 0;
      for (unsigned long oid : clss) {
        unsigned long nid = (clssID << BitsShiftCount) + eindex;
        eindex++;
        NodeID2CNodeID[oid] = nid;
      }
      clssID++;
    } else {
      smallEqualClss.push_back(clss);
//      sum.insert(sum.end(), clss.begin(), clss.end());
    }
  }
//  unsigned int outerindex = 0;
//  for (unsigned long oid : sum) {
//    unsigned long nid = (clssID << 6) + (outerindex);
//    nid2clustering[oid] = nid;
//    outerindex++;
//    if (outerindex == 64) {
//      clssID++;
//      outerindex = 0;
//    }
//  }
  unsigned int outerindex = 0;
  for (std::vector<unsigned long> clss : smallEqualClss) {
    for (unsigned long oid : clss) {
      unsigned long nid = (clssID << BitsShiftCount) + (outerindex);
      NodeID2CNodeID[oid] = nid;
      outerindex++;
      if (outerindex == EquivalenceClassesElementNum) {
        clssID++;
        outerindex = 0;
      }
    }
  }

//  for (auto dd : NodeID2CNodeID) {
//    clusteringFile<< dd.first << ","<< dd.second << std::endl;
//  }
  clusteringFile.close();

  ////////////////////
  if (dir) {
    for (TNGraph::TNodeI NI = Graph1->BegNI(); NI < Graph1->EndNI(); NI++) {
      int nid = NI.GetId();
      std::vector<TEdgeTuple>& AttrFltIntKV =  Graph1()->GetSortedAttrByNode(nid);
      unsigned long size = AttrFltIntKV.size();
      for (int EdgeIndex = 0; EdgeIndex < size; EdgeIndex++) {
        TEdgeTuple& tpl = AttrFltIntKV[EdgeIndex];
        tpl.cto = NodeID2CNodeID[tpl.to];
      }
    }
  }else {
    for (TUNGraph::TNodeI NI = UGraph1->BegNI(); NI < UGraph1->EndNI(); NI++) {
      int nid = NI.GetId();
      std::vector<TEdgeTuple>& AttrFltIntKV =  UGraph1()->GetSortedAttrByNode(nid);
      unsigned long size = AttrFltIntKV.size();
      for (int EdgeIndex = 0; EdgeIndex < size; EdgeIndex++) {
        TEdgeTuple& tpl = AttrFltIntKV[EdgeIndex];
        tpl.cto = NodeID2CNodeID[tpl.to];
      }
    }
  }
  
  
  std::ofstream oFile;
oFile.open((genbase+name+ostr).c_str(), std::ios::out | std::ios::trunc);

oFile<< "NID"<< ","
  << "NumberOfNode" << ","
  << "NumberOfEdge" << ","
 << "SSHeapItem" << ","
 << "DJHeapItem" << ","
 <<  "dj/ss"<< ","
  <<  "SHeapWidth"<< ","
  <<  "DHeapWidth"<< ","
  << "SS-CLOCKS" << "," << "SS-time" << ","
  << "DJ-CLOCKS" << "," << "DJ-time" << std::endl;
 index = 0;
if (dir) {
  for (TNGraph::TNodeI NI = Graph1->BegNI(); NI < Graph1->EndNI(); NI++) {
      if(randomNodeIDIndex[index] == 1){
        int nid = NI.GetId();
//        NodeID2CNodeID.clear();
//      TSnap::nodeID2ClusteringID(NodeID2CNodeID, nid, Graph1());
//        for (auto dd : NodeID2CNodeID) {
//          clusteringFile<< dd.first << ","<< dd.second << std::endl;
//        }
//        clusteringFile.close();
        
          DijkstraStat.restat();
          SmallStepStat.restat();
          
          DijkstraStat.restat();
          SmallStepStat.restat();
          auto graph = Graph1();
        testHeapWithClustering(graph, nid,
               &smallStepTime, &DijkstraTime,
               DijkstraStat, SmallStepStat, NodeID2CNodeID);
        oFile << nid << ","
          << SmallStepStat.getNodeCount() << ","
          << DijkstraStat.getEdgeNum()/SmallStepStat.getNodeCount() << ","
          << (double)SmallStepStat.getHeapItem()/SmallStepStat.getNodeCount() << ","
          << (double)DijkstraStat.getHeapItem()/SmallStepStat.getNodeCount() << ","
          << (double)((double)DijkstraTime/(double)smallStepTime)<< ","
          << SmallStepStat.getUseless()<< ","
        << DijkstraStat.getUseless()<< ","
          << smallStepTime << "," << (double)smallStepTime/CLOCKS_PER_SEC << ","
          << DijkstraTime << "," << (double)DijkstraTime/CLOCKS_PER_SEC
           << std::endl;
        
      }
      index++;
  }
}else {
  for (TUNGraph::TNodeI NI = UGraph1->BegNI(); NI < UGraph1->EndNI(); NI++) {
      if(randomNodeIDIndex[index] == 1){
        int nid = NI.GetId();
//        nid = 234;
//        NodeID2CNodeID.clear();
//      TSnap::nodeID2ClusteringID(NodeID2CNodeID, nid, UGraph1());
//        for (auto dd : NodeID2CNodeID) {
//          clusteringFile<< dd.first << ","<< dd.second << std::endl;
//        }
//        clusteringFile.close();
        
          DijkstraStat.restat();
          SmallStepStat.restat();
          
          DijkstraStat.restat();
          SmallStepStat.restat();
        testHeapWithClustering(UGraph1(), nid,
               &smallStepTime, &DijkstraTime,
               DijkstraStat, SmallStepStat, NodeID2CNodeID);
     oFile << nid << ","
        << SmallStepStat.getNodeCount() << ","
        << DijkstraStat.getEdgeNum()/SmallStepStat.getNodeCount() << ","
          << (double)SmallStepStat.getHeapItem()/SmallStepStat.getNodeCount() << ","
          << (double)DijkstraStat.getHeapItem()/SmallStepStat.getNodeCount()  << ","
        << (double)((double)DijkstraTime/(double)smallStepTime)<< ","
        << SmallStepStat.getUseless()<< ","
        << DijkstraStat.getUseless()<< ","
          << smallStepTime << "," << (double)smallStepTime/CLOCKS_PER_SEC << ","
          << DijkstraTime << "," << (double)DijkstraTime/CLOCKS_PER_SEC
          << std::endl;
      }
      index++;
    
  }
}

oFile.close();
  
}

int main(int argc, char* argv[]) {
  for (int i=0; i < 1; i++) {
     
//  produceClusteringID("bio-CE-GN", false);
//  produceClusteringID("bio-CE-CX", false);//11
//  produceClusteringID("bio-DM-CX", false);//3.6
//  produceClusteringID("bio-HS-CX", false);//3.6
//  produceClusteringID("bio-SC-HT", false);//2.2
//
//  produceClusteringID("bio-human-gene1", false);
////  produceClusteringID("bio-human-gene2", false);
//  produceClusteringID("bio-mouse-gene", false);//5.x
//  produceClusteringID("bio-WormNet-v3", false);//11

//  produceClusteringID("USairport500", true);//11
//  produceClusteringID("OClinks_w_chars", true);//11
//  produceClusteringID("OClinks_w", true);//11
//
//  produceClusteringID("celegans_n306", true);//diff count = 1
//  produceClusteringID("Cross_Parker-Consulting_info", true);//11
//  produceClusteringID("Cross_Parker-Consulting_value", true);//11
//  produceClusteringID("Cross_Parker-Manufacturing_aware", true);//11
//  produceClusteringID("Cross_Parker-Manufacturing_info", true);//11
//  produceClusteringID("Freemans_EIES-1_n48", true);//11
//  produceClusteringID("Freemans_EIES-2_n48", true);//11
  }
//
//  produceClusteringID("datagen-7_5-fb", false);//
//  produceClusteringID("datagen-7_6-fb", false);//
//  produceClusteringID("datagen-7_7-zf", false);//
//  produceClusteringID("datagen-7_8-zf", false);//
//  produceClusteringID("datagen-7_9-fb", false);//
  produceClusteringID("datagen-8_0-fb", false);// 无差别堆操作 不行
  produceClusteringID("datagen-8_1-fb", false);// 无差别堆操作 不行
  
//  produceClusteringID("datagen-8_2-zf", false);//太大了 启用虚拟缓存了
//  produceClusteringID("datagen-9_0-fb", false);//太大了 崩溃
  
  
//    compairSmallStepAndDijkstras("bio-CE-GN", false);
//    compairSmallStepAndDijkstras("bio-CE-CX", false);//11
//    compairSmallStepAndDijkstras("bio-DM-CX", false);//3.6
//    compairSmallStepAndDijkstras("bio-HS-CX", false);//3.6
//    compairSmallStepAndDijkstras("bio-SC-HT", false);//2.2
//
//////    compairSmallStepAndDijkstras("bio-human-gene1", false);
////    compairSmallStepAndDijkstras("bio-human-gene2", false);
//////    compairSmallStepAndDijkstras("bio-mouse-gene", false);//5.x
//    compairSmallStepAndDijkstras("bio-WormNet-v3", false);//11
////
//    compairSmallStepAndDijkstras("USairport500", true);//11
//    compairSmallStepAndDijkstras("OClinks_w_chars", true);//11
//    compairSmallStepAndDijkstras("OClinks_w", true);//11
//
//    compairSmallStepAndDijkstras("celegans_n306", true);//11
//    compairSmallStepAndDijkstras("Cross_Parker-Consulting_info", true);//11
//    compairSmallStepAndDijkstras("Cross_Parker-Consulting_value", true);//11
//    compairSmallStepAndDijkstras("Cross_Parker-Manufacturing_aware", true);//11
//    compairSmallStepAndDijkstras("Cross_Parker-Manufacturing_info", true);//11
//    compairSmallStepAndDijkstras("Freemans_EIES-1_n48", true);//11
//    compairSmallStepAndDijkstras("Freemans_EIES-2_n48", true);//11

//  https://ldbcouncil.org/benchmarks/graphalytics/
//  zstd -d datagen-8_1-fb.tar.zst -o datagen-8_1-fb.tar
  //https://stats.blue/Stats_Suite/multiple_linear_regression_calculator.html
  
//    compairSmallStepAndDijkstras("datagen-7_5-fb", false);//11
//    compairSmallStepAndDijkstras("datagen-7_6-fb", false);//11
//    compairSmallStepAndDijkstras("datagen-7_7-zf", false);//11
//    compairSmallStepAndDijkstras("datagen-7_8-zf", false);//11
//    compairSmallStepAndDijkstras("datagen-7_9-fb", false);//11
//////
////
  


  return 0;
}
