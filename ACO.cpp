//450
#include<bits/stdc++.h>
#include<vector>
#include<map>
using namespace std;
int main(int argc, char* argv[]){
    ios::sync_with_stdio(0);
    cin.tie(0);
    int city;
    int point[100][2];
    char *filename = argv[1];             //儲存城市座標的檔案名稱(ex : .\\dt1\\point.txt)
    ifstream file;
    file.open(filename);                  //開檔
    if(!file.is_open()){                  //檔案無法開啟
        cout << "can not open " << filename << endl;
        return EXIT_FAILURE;
    }
    while(file >> city){                  
        file >> point[city][0] >> point[city][1];        //讀取每個城市的座標
    }
    double distance[city+1][city+1], dis_reciprocal[city+1][city+1], pheromones[city+1][city+1];
    file.close();
    for(int i=1;i<=city;i++){                  //計算每個城市和其他城市的距離
        for(int j=i+1;j<=city;j++){
            distance[i][j] = distance[j][i] = sqrt(pow(point[i][0] - point[j][0],2) + pow(point[i][1] - point[j][1],2));
            dis_reciprocal[i][j] = dis_reciprocal[j][i] = (1/(distance[i][j]));             //距離倒數
            //pheromones[i][j] = pheromones[j][i] = 1;                                      //初始費洛蒙濃度為1
        }
    }
    int s[city];           //儲存當前走訪順序
    int min_index;
    double ans = 0;
    int cur_dis = 0;                 //兩點之間的距離
    int pre_index = 1;               //前一個點的編號
    double min_dis = 0;                 //前往下一個點的最小距離
    double total_dis = 0;               //總共經過的距離
    bool v[city+1]={0};      //紀錄是否已走過該城市
    v[1] = 1;                //第一個點紀錄為已走過
    int seq1[city]={1};          //紀錄走訪順序
    for(int i=1;i<city;i++){     //尋找下一個最近的點
        min_dis = 0;
        for(int j=1;j<=city;j++){
            if(!v[j]){
                cur_dis = distance[pre_index][j];
                if(!min_dis || (cur_dis < min_dis)){              //紀錄最短距離以及編號
                    min_dis = cur_dis;
                    min_index = j;
                }
            }
        }
        v[min_index] = 1;                      //將找到的點設定為已走過
        total_dis += min_dis;                      //加上該段距離
        seq1[i] = min_index;                       //將找到的編號填入
        pre_index = min_index;                     
    }
    total_dis += distance[pre_index][1];    //加上最後一個點到起點的距離
    if(!ans || total_dis < ans){          //判斷是否小於原本找到的最短距離
        for(int k=0;k<city;k++){
            s[k] = seq1[k];
        }
        ans = total_dis;
    }
    double nn = ans;
    double nn_dis = (1/ans) / city;
    double p = 0.1, a=0.1, b=2, q0 = 0.9, pr, psum, acc, rd, delta, min, aver, argmax,local_min,best_ans=10000;
    int ant_num = 20, max = 0, argmax_index;
    vector<double> prob;
    vector<int> pos;
    vector<int> ans_seq;
    vector<int> best;
    bool foo;
    double dis[ant_num];
    vector<int> seq[ant_num];
    int visit[city+1];
    srand(time(NULL));
    for(int i=0;i<atoi(argv[3]);i++){               // n run
        ans = 1000000000;
        for(int j=1;j<city;j++){
            for(int k=j+1;k<=city;k++){
                pheromones[j][k] = pheromones[k][j] = nn_dis;                                      //初始費洛蒙濃度
                //cout << pheromones[j][k] << endl;
                //cout << pheromones[j][k] << endl;
            }
        }
        // for(int j=0;j<city-1;j++){
        //     pheromones[s[j]][s[j+1]] *= 2;
        // }
        // pheromones[s[city-1]][s[0]] *= 2;
        for(int j=0;j<atoi(argv[2]);j++){           // n iteations in 1 run
            /*if(j%100 == 0)
                cout << j << endl;*/
            local_min = 10000000;
            for(int k=0;k<ant_num;k++){           // n ants in 1 iteration
                
                dis[k] = 0;                        //總距離
                seq[k].push_back(rand() % city + 1);                     //隨機起點  
                //seq[k] = rand() % city + 1;
                memset(visit,0,sizeof(visit));   
                visit[seq[k][0]] = 1;                    
                for(int x=1;x<city;x++){
                    rd = double(rand()) / double(RAND_MAX);                //0-1的隨機數(決定探索路徑的尋找方法)
                    argmax=-1;
                    if(rd > q0){
                        psum = 0;
                        for(int l=1;l<=city;l++){
                            if(!visit[l]){              //到編號l城市的可能性
                                pos.push_back(l);
                                pr = pheromones[seq[k][x-1]][l] * pow(dis_reciprocal[seq[k][x-1]][l], b);
                                //cout << l << " " << pow(dis_reciprocal[seq[k][x-1]][l], b) << endl;
                                //cout << pr << endl;
                                psum += pr;
                                prob.push_back(pr);
                            }    
                        }
                        rd = double(rand()) / double(RAND_MAX);                //0-1的隨機數(決定要走哪條路)
                        acc = 0;
                        for(int l=0;l<prob.size();l++){
                            prob[l] /= psum;                                 //標準化成0-1的區間
                            acc += prob[l];
                            if(acc >= rd){                                   //總和 > 隨機數
                                dis[k] += distance[seq[k][x-1]][pos[l]];     //將距離加上
                                seq[k].push_back(pos[l]);                          //走到編號l的城市
                                visit[pos[l]] = 1;                           //設為已走過
                                break;
                            }
                            if(l == (prob.size()-1)){
                                dis[k] += distance[seq[k][x-1]][pos[l]];     //將距離加上
                                seq[k].push_back(pos[l]);                          //走到編號l的城市
                                visit[pos[l]] = 1;
                                break;
                            }
                        }
                        pos.clear();
                        prob.clear();
                    }
                    
                    else{
                        for(int l=1;l<=city;l++){
                            if(!visit[l]){              //到編號l城市的可能性
                                pr = pheromones[seq[k][x-1]][l] * pow(dis_reciprocal[seq[k][x-1]][l], b);
                                /*if(pr == 0){
                                    cout << pheromones[seq[k][x-1]][l] << endl;
                                    cout << pow(dis_reciprocal[seq[k][x-1]][l], b) << endl;
                                }*/
                                //cout << pheromones[seq[k][x-1]][l] / pow(dis_reciprocal[seq[k][x-1]][l], b) << endl;
                                if(argmax < pr){
                                    argmax = pr;
                                    argmax_index = l;
                                }
                            }    
                        }
                        seq[k].push_back(argmax_index);
                        dis[k] += distance[argmax_index][seq[k][x-1]];
                        visit[argmax_index] = 1;                 //紀錄走過的點
                    }
                } 
                dis[k] += distance[seq[k][0]][seq[k][city-1]];
                
                if(dis[k] < ans){                               //只有得到當前最佳解的螞蟻可以更新費洛蒙
                    delta = 1 / dis[k];
                    for(int n=1;n<city;n++){
                        for(int m=n+1;m<=city;m++){
                            pheromones[n][m] *= (1-a);
                            pheromones[m][n] *= (1-a);
                        }
                    }
                    // //cout << pheromones[seq[k][city-1]][seq[k][0]] << endl;
                    for(int n=0;n<city-1;n++){
                        pheromones[seq[k][n]][seq[k][n+1]] += a*delta;
                        pheromones[seq[k][n+1]][seq[k][n]] += a*delta;
                    }
                    //cout << "b  " << pheromones[seq[k][city-1]][seq[k][0]] << endl;
                    pheromones[seq[k][city-1]][seq[k][0]] += a*delta;
                    pheromones[seq[k][0]][seq[k][city-1]] += a*delta;
                    //cout << pheromones[seq[k][city-1]][seq[k][0]] << endl;
                    //cout << ans << endl;
                    //cout << "a  " <<  pheromones[seq[k][city-1]][seq[k][0]] << endl;
                    ans = dis[k];
                    ans_seq = seq[k];
                }
                
                pheromones[seq[k][1]][seq[k][0]] = pheromones[seq[k][1]][seq[k][0]]*(1-p) + p*nn_dis;       //本地更新
                pheromones[seq[k][0]][seq[k][1]] = pheromones[seq[k][0]][seq[k][1]]*(1-p) + p*nn_dis;
                
            }
            
            for(int k=0;k<ant_num;k++){         
                    seq[k].clear();
            }
            
        }
        
        
        cout << ans << endl;
        if(ans < best_ans){
            best_ans=ans;
            best = ans_seq;
        }
        aver += ans;
    }
    cout << "average : " << aver / atoi(argv[3]) << endl;
    cout << "best : " << best_ans << endl;

    ofstream file1;
    string str(argv[1]);
    file1.open(str.substr(str.find("\\") + 2,str.find("\\point")-str.find("\\") - 3) + "\\ans_" + str.substr(str.find("\\") + 2,str.find("\\point")-str.find("\\") - 3) + ".txt");  //創建新的txt檔
    file1 << "distance: " << best_ans << endl;
    for(int i=0;i<city;i++){
        file1 << best[i] << endl;
    }
    file1.close();
}