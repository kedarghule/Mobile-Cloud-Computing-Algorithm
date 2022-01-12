#include <iostream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <ctime>
using namespace std;

struct Node
{
	int node_no;
	vector<Node*> predecessor;
	vector<Node*> successor;
	int T_l_min;
	int T_re;
	char execution;
	double weight;
	double avg_T_l_min;
	double priority;
	char scheduled_is;
	int core;
	double RT;
	double FT;
	double FT_cloud;
	double ST;
	vector<char> TM_execution;
	vector<double> TM_FT;
	vector<double> TM_ST;
	vector<double> TM_core;
	vector<double> TM_FT_cloud;
	vector<double> TM_RT;
};

struct Result_Sets
{
	vector<vector<int> > results;
};

class TaskGraph
{
	Node* initialize_Node(int vertex_val)
	{
		Node* new_node = new Node;
		new_node->node_no = vertex_val;
		new_node->predecessor = vector<Node*>();
		new_node->successor = vector<Node*>();
		new_node->T_l_min = 0;
		new_node->T_re = 0;
		new_node->execution = 'L';
		new_node->weight = 0.0;
		new_node->avg_T_l_min=0.0;
		new_node->priority=0.0;
		new_node->scheduled_is = 'N';
		new_node->RT = 0;
		new_node->FT = 0;
		new_node->FT_cloud = 0;
		new_node->ST = 0;
		new_node->core = 0;
		return new_node;
	}

    	int N_Tasks;
    	int N_Cores;
    	int T_max = 27; // Given Constraint
    	int N_local = 0; // Number of tasks on local cores
    	double E_initial_schedule = 0.0;
    	double T_initial_schedule = 0.0;	
    	
	public:
	    Node **head;                
	    vector<Node*> Task_Node;
	    vector<Node*> descendingOrderedTasks;
	    vector<Node*> LocalTasks;
	    vector<double> Core_Time_Pointer;
	    vector<double> Cloud_Time_Pointer;
	    vector<vector<int> > Set_k;
	    vector<double> E_total;
	    vector<double> T_total;
	    vector<Result_Sets*> Set_Results_K;
	    vector<vector<int> > MASTER_Results;
	    int MASTER_Choice_Index;
	    int MASTER_Energy_Index;
	    
	    // Constructor to initialize the Graph of Tasks
	    template<typename T>
	    TaskGraph(vector<T> &graph_edges, int n_edges, int n_tasks, int n_cores)  
	    {
	    	N_Tasks = n_tasks;
	    	N_Cores = n_cores;
	    	Core_Time_Pointer.resize(n_cores, 0);
	    	Cloud_Time_Pointer.resize(3, 0);
	    	Set_k.resize(n_cores+1);
	    	
	    	for(int i = 0; i<n_tasks; i++)
	    	{
	    		Node* new_node = new Node;
	    		new_node = initialize_Node(i+1);
	    		Task_Node.push_back(new_node);
	    	}
	    	
	    	for(int i=0;i<n_tasks;i++)
	    	{
	    		for(int j=0;j<n_edges;j++)
		    	{
		    		int start_vertex = graph_edges.at(j).first;
		    		if(Task_Node[i]->node_no == start_vertex)
		    		{
		    			int end_vertex = graph_edges.at(j).second;
		    			Task_Node[i]->successor.push_back(Task_Node[end_vertex-1]);
		    			Task_Node[end_vertex-1]->predecessor.push_back(Task_Node[start_vertex-1]);
		    		}
			}
	    	}
	    }
	    
	    // Function to display initial graph details
	    void displayGraphDetails()
	    {
	    	for(int i = 0; i<N_Tasks; i++)
	    	{
	    		cout<<"Node "<<i+1<<" has task number "<<Task_Node[i]->node_no<<endl;
	    		cout<<"Its predecessors are :"<<endl;
	    		for(int j=0;j<Task_Node[i]->predecessor.size();j++)
	    			cout<<Task_Node[i]->predecessor[j]->node_no<<" ";
			cout<<endl;
			cout<<"Its successors are :"<<endl;
			for(int j=0;j<Task_Node[i]->successor.size();j++)
				cout<<Task_Node[i]->successor[j]->node_no<<" ";
			cout<<endl;
					
			cout<<"---------------------------------"<<endl;
	    	}
	    }
	    
	    // STEP 1:PHASE 1 -Function for primary assignment
	    void primary_assignment(int *T_local, int T_cloud[], int n_tasks, int n_cores)
	    {
	    	for(int i = 0; i<n_tasks; i++)
	    	{
	    		int sum = 0;
   	  		Task_Node[i]->T_re = accumulate(T_cloud, T_cloud+3, sum);
   	  		int minimum = *((T_local+i*n_cores) + 0);
			sum =0;
   	  		for(int j=0;j<n_cores;j++)
   	  		{
   	  			if(minimum>*((T_local+i*n_cores) + j)) 
   	  			{
					minimum=*((T_local+i*n_cores) + j);
			      	}
			      	sum+=*((T_local+i*n_cores) + j);
   	  		}
   	  		Task_Node[i]->T_l_min = minimum;
   	  		Task_Node[i]->avg_T_l_min = (double)sum/n_cores;
   	  		   	  		
   	  		if(Task_Node[i]->T_re<Task_Node[i]->T_l_min)
   	  			Task_Node[i]->execution = 'C';
   	  		else
   	  			Task_Node[i]->execution = 'L';
	    	}
	    }
	    
	    // Function to get the weight of a task
	    double get_weight(int i)
	    {
	    	if(Task_Node[i]->execution=='C')
	    		return Task_Node[i]->T_re;
	    	else
	    		return Task_Node[i]->avg_T_l_min;
	    }
	    
	    // Function to get the maximum priority among the successors of a task
	    double get_max_priority(int i)
	    {
	    	double max = Task_Node[i]->successor[0]->priority;
	    	for(int j=1;j<Task_Node[i]->successor.size();j++)
	    	{
	    		if(max<Task_Node[i]->successor[j]->priority)
	    			max = Task_Node[i]->successor[j]->priority;
	    	}
	    	return max;
	    }
	    
	    // STEP 1:PHASE 2 -Recursive function to find Task Priorities
	    void task_prioritizing(int i) // i is initially passed as n_tasks
	    {
	    	if(Task_Node[i]->node_no == 1)
	    	{
	    		Task_Node[i]->weight=get_weight(i);
	    		Task_Node[i]->priority = (double)Task_Node[i]->weight + get_max_priority(i);
	    		return;
	    	}
	    	else
	    	{
	    		if(Task_Node[i]->successor.empty())
	    		{
	    			Task_Node[i]->weight=get_weight(i);
	    			Task_Node[i]->priority=Task_Node[i]->weight;
	    		}
	    		else
	    		{
	    			Task_Node[i]->weight=get_weight(i);
	    			Task_Node[i]->priority = (double)Task_Node[i]->weight + get_max_priority(i);
	    		}
	    	}
	    	task_prioritizing(i-1);
	    }
	    
	    // Check if the task's predecessors are scheduled or not
	    bool getPredScheduledValue(int i)
	    {
	    	bool checkPredScheduled = true;
	    	for(int j=0;j<Task_Node[i]->predecessor.size();j++)
	    		if(Task_Node[i]->predecessor[j]->scheduled_is == 'N')
	    			checkPredScheduled = false;
	    	return checkPredScheduled;
	    }
	    
	    // STEP 1:PHASE 3 - Function for Execution Unit Selection
	    void execution_unit_selection(int *T_local, int T_cloud[], float Pk[], float Ps)
	    {
	    	// Sorting tasks in descending order of priorities
	    	descendingOrderedTasks = Task_Node;
	    	sort(descendingOrderedTasks.begin(), descendingOrderedTasks.end(), [](Node* a, Node* b) {
			return a->priority > b->priority;
		});
		
		for(int i=0;i<N_Tasks;i++)
		{
			int node_number = descendingOrderedTasks[i]->node_no;
			bool checkPredScheduled = getPredScheduledValue(node_number-1);
			if(checkPredScheduled && descendingOrderedTasks[i]->scheduled_is == 'N')
			{
				double sum_cloud_wr = 0.0;
				double max_FT_c_j = 0.0;
				double max_FT_j = 0.0;
				// Check if tasks has predecessors or not and get maximum finish time values from the predecessors
				if (descendingOrderedTasks[i]->predecessor.size() != 0)
				{
					max_FT_c_j = descendingOrderedTasks[i]->predecessor[0]->FT_cloud;
					max_FT_j = descendingOrderedTasks[i]->predecessor[0]->FT;
					for(int j=1;j<descendingOrderedTasks[i]->predecessor.size();j++)
					{
						if(max_FT_c_j<descendingOrderedTasks[i]->predecessor[j]->FT_cloud)
							max_FT_c_j=descendingOrderedTasks[i]->predecessor[j]->FT_cloud;
						if(max_FT_j<descendingOrderedTasks[i]->predecessor[j]->FT)
							max_FT_j=descendingOrderedTasks[i]->predecessor[j]->FT;
					}
				}
				// Calculating time taken by task if executed on cloud
				double FT_ws_i = max(Cloud_Time_Pointer[0],max_FT_j) + T_cloud[0];
				double RT_c_i = max(FT_ws_i, max_FT_c_j); // equation 5
				double FT_c_i = RT_c_i + T_cloud[1];
				sum_cloud_wr = FT_c_i + T_cloud[2];
				
				 // Task Assigned to Local during Primary Assignment
				if(descendingOrderedTasks[i]->execution == 'L')
				{
					vector<double> core_sum(N_Cores,0);
					for(int j=0;j<N_Cores;j++)
					{
						core_sum[j] = max(Core_Time_Pointer[j], max_FT_j) + *((T_local+(node_number-1)*N_Cores) + j);
					}
					auto minimum_core_time_iterator = min_element(begin(core_sum), end(core_sum));
					int core = distance(begin(core_sum), minimum_core_time_iterator);
					double minimum_core_time = core_sum[core];
					
					if(sum_cloud_wr<minimum_core_time)
					{
						//Task is scheduled on cloud and values are updated
						descendingOrderedTasks[i]->execution='C';
						descendingOrderedTasks[i]->scheduled_is='Y';
						descendingOrderedTasks[i]->FT=sum_cloud_wr;
						descendingOrderedTasks[i]->RT=max_FT_j;
						descendingOrderedTasks[i]->core=-1;
						descendingOrderedTasks[i]->FT_cloud = FT_c_i;
						Cloud_Time_Pointer[0] = FT_ws_i; // Ready time for wireless send for next task
						E_initial_schedule = E_initial_schedule + T_cloud[0]*Ps;
						Set_k[0].push_back(node_number);
						int T_wr_sum = 0;
   	  					descendingOrderedTasks[i]->ST = descendingOrderedTasks[i]->FT - accumulate(T_cloud, T_cloud+3, T_wr_sum);
					}
					else
					{
						//Task is scheduled on a core and values are updated
						descendingOrderedTasks[i]->execution='L';
						descendingOrderedTasks[i]->scheduled_is='Y';
						descendingOrderedTasks[i]->FT=minimum_core_time;
						descendingOrderedTasks[i]->RT=max_FT_j;
						descendingOrderedTasks[i]->core=core+1;
						Core_Time_Pointer[core] = minimum_core_time;
						E_initial_schedule = E_initial_schedule + *((T_local+(node_number-1)*N_Cores) + core)*Pk[core];
						N_local++;
						Set_k[core+1].push_back(node_number);
						descendingOrderedTasks[i]->ST = descendingOrderedTasks[i]->FT - *((T_local+(node_number-1)*N_Cores) + core);
						LocalTasks.push_back(descendingOrderedTasks[i]);
					}
				}
				else // Task Assigned to Cloud during Primary Assignment
				{
					//Task is scheduled on cloud and values are updated
					descendingOrderedTasks[i]->execution='C';
					descendingOrderedTasks[i]->scheduled_is='Y';
					descendingOrderedTasks[i]->FT=sum_cloud_wr;
					descendingOrderedTasks[i]->RT=max_FT_j;
					descendingOrderedTasks[i]->core=-1;
					descendingOrderedTasks[i]->FT_cloud = FT_c_i;
					Cloud_Time_Pointer[0] = FT_ws_i; // Ready time for wireless send for next task
					E_initial_schedule = E_initial_schedule + T_cloud[0]*Ps;
					Set_k[0].push_back(node_number);
					int T_wr_sum = 0;
   	  				descendingOrderedTasks[i]->ST = descendingOrderedTasks[i]->FT - accumulate(T_cloud, T_cloud+3, T_wr_sum);
				}	
			}
		}
		
		// Displaying the results of Initial Scheduling
		cout<<"************ INITIAL SCHEDULING RESULTS **************"<<endl;
		for(int i = 0; i<N_Tasks; i++)
	    	{
	    		cout<<"TASK NUMBER = "<<Task_Node[i]->node_no<<endl;
	    		cout<<"Its predecessors are : ";
	    		for(int j=0;j<Task_Node[i]->predecessor.size();j++)
	    			cout<<Task_Node[i]->predecessor[j]->node_no<<" ";
			cout<<endl;
			cout<<"Its successors are : ";
			for(int j=0;j<Task_Node[i]->successor.size();j++)
				cout<<Task_Node[i]->successor[j]->node_no<<" ";
			cout<<endl;
			cout<<"Task is assigned to :"<<Task_Node[i]->execution<<endl;
			cout<<"Priority Value :"<<Task_Node[i]->priority<<endl;
			cout<<"Scheduling Status :"<<Task_Node[i]->scheduled_is<<endl;
			cout<<"Core :"<<Task_Node[i]->core<<endl;
			cout<<"Ready Time :"<<Task_Node[i]->RT<<endl;
			cout<<"Start Time :"<<Task_Node[i]->ST<<endl;
			cout<<"Finish Time :"<<Task_Node[i]->FT<<endl;			
			cout<<"---------------------------------"<<endl;
			if(descendingOrderedTasks[i]->successor.size() ==0)
				if(T_initial_schedule<descendingOrderedTasks[i]->FT)
					T_initial_schedule = descendingOrderedTasks[i]->FT;
							
	    	}
		
	    	cout<<"Energy Consumed : "<<E_initial_schedule<<endl;
	    	cout<<"Application Completion Time : "<<T_initial_schedule<<endl;    	
	    	cout<<"The initial scheduling order is : "<<endl;
	    	for(int i = 0;i<Set_k.size();i++)
	    	{
	    		if(i==0)
	    			cout<<"On cloud, we have following tasks : ";
	    		else
	    			cout<<"In kernel "<<i<<", we have following tasks : ";
	    		for(int j = 0;j<Set_k[i].size();j++)
	    			cout<<Set_k[i][j]<<" ";
	    		cout<<endl;
	    			
	    	}
	    	cout<<"*********************************************"<<endl;
	    }
	    
	    // STEP 2 - Function for Task Migration
	    void Task_Migration(int *T_local, int T_cloud[], float Pk[], float Ps)
	    {
		MASTER_Energy_Index = 0;
	      	int TM_counter=0;
	    	
	    	int while_loop_counter = 0;
	    	
	    	// A While Loop that runs your Task Migration Algorithm to optimize results
	    	while(true)
	    	{
		    	while_loop_counter++;
		    	//Outer Loop
		    	for(int i = 0; i<N_local; i++)
		    	{
		    		vector<int> K(N_Cores+1) ; // vector with 4 ints.
				iota (begin(K), end(K), 0); // Fill with 0, 1, ..., 99.
				
				remove(K.begin(),K.end(),LocalTasks[i]->core);
				K.resize(N_Cores);
							
				for(int j=0;j<N_Cores;j++)
				{
					// TASK MIGRATION ALGORITHM
					
					// Reset/Set the Cloud and Core Pointer vectors
					vector<double> TM_Cloud_Time_Pointer;
		    			vector<double> TM_Core_Time_Pointer;
					TM_Core_Time_Pointer.resize(N_Cores, 0);
		    			TM_Cloud_Time_Pointer.resize(3, 0);
		    			
		    			// Create Set K New
		    			vector<vector<int> >  Set_k_new;
					Set_k_new.resize(N_Cores+1);
					for(int p=0;p<N_Cores+1;p++)
					{
						Set_k_new[p].resize(Set_k[p].size());
					}
					
					for(int p=0;p<Set_k_new.size();p++)
					{
						for(int q=0;q<Set_k_new[p].size();q++)
							Set_k_new[p][q] = Set_k[p][q];
					}
								
					// Remove current ith task in LocalTasks from the set
					remove(Set_k_new[LocalTasks[i]->core].begin(),Set_k_new[LocalTasks[i]->core].end(),LocalTasks[i]->node_no);
					Set_k_new[LocalTasks[i]->core].resize(Set_k_new[LocalTasks[i]->core].size()-1);
					
					// Assign the removed task a new position
					if(Set_k_new[K[j]].size() == 0)
					{
						Set_k_new[K[j]].push_back(LocalTasks[i]->node_no);
					}
					else
					{
						int pos = 0;
						for(int p=0;p<Set_k_new[K[j]].size();p++)
						{
							if(LocalTasks[i]->RT>Task_Node[Set_k_new[K[j]][p]-1]->ST)
								pos = p+1;
							else
							{
								pos=p;
								break;
							}
						}
						Set_k_new[K[j]].insert(Set_k_new[K[j]].begin()+pos, LocalTasks[i]->node_no);
					}
					
				    	// Add new scheduling set in results set
				    	Result_Sets* new_node = new Result_Sets;
		    			new_node->results = Set_k_new;
		    			Set_Results_K.push_back(new_node);
				    	
				    	// Update the cores of the all the tasks to the new tasks cores in the above permutation
				    	for(int q = 0;q<Set_k_new.size();q++)
				    		for(int r=0;r<Set_k_new[q].size();r++)
				    		{
				    			int task_number = Set_k_new[q][r];
				    			Task_Node[task_number-1]->TM_core.push_back(q);
				    		}
				    	
				    	// Initialize vector ready1
					vector<int> ready1;
					ready1.resize(N_Tasks);
					for(int q=0;q<ready1.size();q++)
					{
						Task_Node[q]->scheduled_is = 'N';
						if(getPredScheduledValue(q) == true)
							ready1[q]=0;
						else
						{
							int count=0;
							for(int r=0;r<Task_Node[q]->predecessor.size();r++)
							{
								if(Task_Node[q]->predecessor[r]->scheduled_is == 'N')
									count++;
							}
							ready1[q]=count;
						}
					}
					
					// Initialize vector ready2
					vector<int> ready2;
					ready2.resize(N_Tasks,1);
					ready2[0] = 0;
					for(int p = 0;p<Set_k_new.size();p++)
				    	{
				    		int task_number;
				    		if(Set_k_new[p].size()!=0)
				    		{
				    			task_number = Set_k_new[p][0];
				    			ready2[task_number-1] = 0;
				    		}
				    	}
					
					// Initialize the LIFO Stack
					vector<int> Stack_schedule;
					for(int p=0;p<ready2.size();p++)
					{
						if(ready1[p] == 0 && ready2[p]==0)
							Stack_schedule.push_back(Task_Node[p]->node_no);
					}
					
					double E_TM = 0.0;
					double T_TM = 0.0;
					
					// Scheduling of tasks. Loop runs until stack gets empty
					while(Stack_schedule.size()!=0)
					{
						int v_i = Stack_schedule.back();
						v_i--;
						Stack_schedule.pop_back();
						
						double sum_cloud_wr = 0.0;
						double max_FT_c_j = 0.0;
						double max_FT_j = 0.0;
						
						if (Task_Node[v_i]->predecessor.size() != 0)
						{
							max_FT_c_j = Task_Node[v_i]->predecessor[0]->TM_FT_cloud[TM_counter];
							max_FT_j = Task_Node[v_i]->predecessor[0]->TM_FT[TM_counter];
							
							for(int p=1;p<Task_Node[v_i]->predecessor.size();p++)
							{
								if(max_FT_c_j<Task_Node[v_i]->predecessor[p]->TM_FT_cloud[TM_counter])
									max_FT_c_j=Task_Node[v_i]->predecessor[p]->TM_FT_cloud[TM_counter];
								if(max_FT_j<Task_Node[v_i]->predecessor[p]->TM_FT[TM_counter])
									max_FT_j=Task_Node[v_i]->predecessor[p]->TM_FT[TM_counter];
							}
							
						}
						
						double FT_ws_i = max(TM_Cloud_Time_Pointer[0],max_FT_j) + T_cloud[0];
						double RT_c_i = max(FT_ws_i, max_FT_c_j); // equation 5
						double FT_c_i = RT_c_i + T_cloud[1];
						sum_cloud_wr = FT_c_i + T_cloud[2];
						int node_number = Task_Node[v_i]->node_no;	
						
						// If Set Assigns task to cloud
						if(Task_Node[v_i]->TM_core[TM_counter] == 0)
						{
							Task_Node[v_i]->TM_execution.push_back('C');
							Task_Node[v_i]->scheduled_is='Y';
							Task_Node[v_i]->TM_FT.push_back(sum_cloud_wr);
							int T_wr_sum = 0;
							Task_Node[v_i]->TM_ST.push_back(Task_Node[v_i]->TM_FT[TM_counter] - accumulate(T_cloud, T_cloud+3, T_wr_sum));
							Task_Node[v_i]->TM_FT_cloud.push_back(FT_c_i);
							TM_Cloud_Time_Pointer[0] = FT_ws_i; // Ready time for wireless send for next task
							E_TM = E_TM + T_cloud[0]*Ps;
							
							Task_Node[v_i]->TM_RT.push_back(max_FT_j);
						}
						// Else set assigns task to core
						else
						{
							int current_core = Task_Node[v_i]->TM_core[TM_counter];
							double minimum_core_time = max(TM_Core_Time_Pointer[current_core-1], max_FT_j)+*((T_local+(node_number-1)*N_Cores) + (current_core-1));
							Task_Node[v_i]->TM_execution.push_back('L');
							Task_Node[v_i]->scheduled_is='Y';
							Task_Node[v_i]->TM_FT.push_back(minimum_core_time);
							
							Task_Node[v_i]->TM_ST.push_back(Task_Node[v_i]->TM_FT[TM_counter] - *((T_local+(node_number-1)*N_Cores) + (current_core-1)));
							Task_Node[v_i]->TM_FT_cloud.push_back(0);
							TM_Core_Time_Pointer[current_core-1] = minimum_core_time;
							E_TM = E_TM + *((T_local+(node_number-1)*N_Cores) + (current_core-1))*Pk[current_core-1];
							Task_Node[v_i]->TM_RT.push_back(max_FT_j);
						}
						int current_core = Task_Node[v_i]->TM_core[TM_counter];
						
						// Update Ready1 and Ready2 for that task to be -1 so that task never gets selected again for schedule
						ready2[v_i] = -1;
						ready1[v_i] = -1;
						
						// Update ready1 for all tasks
						for(int p = 0;p<Task_Node[v_i]->successor.size();p++)
						{
							ready1[Task_Node[v_i]->successor[p]->node_no-1]--;
						}
						
						// Update ready2 for all tasks
						if(Task_Node[v_i]->node_no != Set_k_new[current_core].back())
						{
							for(int p=0;p<Set_k_new[current_core].size();p++)
							{
								if(Set_k_new[current_core][p] == Task_Node[v_i]->node_no)
								{
									int next_local = Set_k_new[current_core][p+1];
									next_local--;
									ready2[next_local] = 0;
								}
							}
						}
						
						// Push in Stack
						for(int p=0;p<ready2.size();p++)
						{
							if(ready1[p] == 0 && ready2[p]==0 && count(Stack_schedule.begin(), Stack_schedule.end(), Task_Node[p]->node_no)==0)
								Stack_schedule.push_back(Task_Node[p]->node_no);
						}					
					}
					
					// Find total time taken
					double T_task_migration_initial = 0.0;
					for(int p=0;p<N_Tasks;p++)
						if(Task_Node[p]->successor.size() ==0)
							if(T_task_migration_initial<Task_Node[p]->TM_FT[TM_counter])
								T_task_migration_initial = Task_Node[p]->TM_FT[TM_counter];
					
					// Push energy and time in a vector of results
					E_total.push_back(E_TM);
					T_total.push_back(T_task_migration_initial);
					TM_counter++;
				}
		    	}

		    	// Make choice for finding the optimal task schedule
		    	double E_choice1, T_choice1;
		    	double E_choice2, T_choice2;
		    	int choice_index1, choice_index2, choice_index; 
		    	int choice1_flag=0; int choice2_flag=0;
		    	double E_red=0.0;
		    	double max_ratio =0.0;
		    	double E_choice, T_choice;
		    	for(int i = MASTER_Energy_Index;i<E_total.size();i++)
		    	{
		    		if(T_total[i]<=T_max && E_total[i]<E_initial_schedule)
		    		{
		    			double ratio = (E_initial_schedule - E_total[i])/(T_total[i] - T_initial_schedule);
		    			if(T_total[i]<=T_initial_schedule && E_red<=E_initial_schedule - E_total[i])
		    			{
			    			E_red = E_initial_schedule - E_total[i];
						E_choice1 = E_total[i];
						T_choice1 = T_total[i];
						choice_index1 = i;
						choice1_flag++;
					}
					else if(max_ratio<ratio)
					{
						max_ratio = ratio;
						E_choice2 = E_total[i];
						T_choice2= T_total[i];
						choice_index2 = i;
						choice2_flag++;
					}					
		    		}
		    		
		    	}
		    	if(choice1_flag)
		    	{
		    		E_choice = E_choice1;
		    		T_choice = T_choice1;
		    		choice_index = choice_index1;
		    	}
		    	else if(choice2_flag)
		    	{
		    		E_choice = E_choice2;
		    		T_choice = T_choice2;
		    		choice_index = choice_index2;
		    	}
		    	
		    	// If energy has not decreased, we break out of the outer while loop
			if(choice1_flag==0 && choice2_flag==0)
			{
				break;
			}
			// else we update the initial values with the optimal values we just got 
			// and run another iteration of the while loop
			else
			{
				MASTER_Results = Set_Results_K[choice_index]->results;
				MASTER_Choice_Index = choice_index;
				Set_k = MASTER_Results;
				E_initial_schedule = E_total[MASTER_Choice_Index];
				T_initial_schedule = T_total[MASTER_Choice_Index];
				MASTER_Energy_Index = E_total.size()-1;
				LocalTasks.clear();
				for(int i=0;i<N_Tasks;i++)
				{
					if(Task_Node[i]->TM_execution[MASTER_Choice_Index] == 'L')
					{
						Task_Node[i]->core = Task_Node[i]->TM_core[MASTER_Choice_Index];
						Task_Node[i]->ST = Task_Node[i]->TM_ST[MASTER_Choice_Index];
						Task_Node[i]->RT = Task_Node[i]->TM_RT[MASTER_Choice_Index];
						LocalTasks.push_back(Task_Node[i]);
					}
				}
				N_local = LocalTasks.size(); 
			}
	
		} // while loop ends here
		
		// Displaying results of Task migration
		cout<<"*************** TASK MIGRATION RESULTS *******************"<<endl;
		cout<<"Scheduling Results for each task after Task Migration are : "<<endl;
		for(int i = 0; i<N_Tasks; i++)
	    	{
	    		cout<<"TASK NUMBER =  "<<Task_Node[i]->node_no<<endl;
	    		cout<<"Its predecessors are : ";
	    		for(int j=0;j<Task_Node[i]->predecessor.size();j++)
	    			cout<<Task_Node[i]->predecessor[j]->node_no<<" ";
			cout<<endl;
			cout<<"Its successors are : ";
			for(int j=0;j<Task_Node[i]->successor.size();j++)
				cout<<Task_Node[i]->successor[j]->node_no<<" ";
			cout<<endl;
			cout<<"Task is assigned to : "<<Task_Node[i]->TM_execution[MASTER_Choice_Index]<<endl;
			cout<<"Core :"<<Task_Node[i]->TM_core[MASTER_Choice_Index]<<endl;
			cout<<"Scheduling Status :"<<Task_Node[i]->scheduled_is<<endl;
			cout<<"Ready Time :"<<Task_Node[i]->TM_RT[MASTER_Choice_Index]<<endl;
			cout<<"Start Time :"<<Task_Node[i]->TM_ST[MASTER_Choice_Index]<<endl;
			cout<<"Finish Time :"<<Task_Node[i]->TM_FT[MASTER_Choice_Index]<<endl;
			cout<<"---------------------------------"<<endl;							
	    	}
		cout<<"Optimal Energy Consumed: "<<E_total[MASTER_Choice_Index]<<endl;
		cout<<"Optimal Application Completion Time : "<<T_total[MASTER_Choice_Index]<<endl;
		cout<<"Scheduling After Task Migration looks like : "<<endl;
		for(int i = 0;i<Set_k.size();i++)
		{
			if(i==0)
	    			cout<<"On cloud, we have following tasks : ";
	    		else
	    			cout<<"In kernel "<<i<<", we have following tasks : ";
		    	for(int j = 0;j<Set_k[i].size();j++)
				cout<<Set_k[i][j]<<" ";
		    	cout<<endl;
		}
		cout<<"****************************************************************"<<endl;
	    }
};

template<typename T>
void printVectorElements(vector<T> &vec)
{
    for (int i = 0; i < vec.size(); i++) 
    {
        cout << "(" << vec.at(i).first << ","<< vec.at(i).second << ")" << "; ";
    }
    cout << endl;
}

int main()
{
	int T_cloud[3] = {3, 1, 1}; // An array containing Ts, Tc and Ti values
	
	int t=10, c=3, edge_no=15; // t denoted number of tasks, c denotes number of core and edge_no denotes number of edges in the graph
	// Values for local cores
	int T_local[t][c] = {{9,7,5}, {8,6,5},{6,5,4},{7,5,3},{5,4,2},{7,6,4},{8,5,3},{6,4,2},{5,3,2},{7,4,2}};
	// Edges in the graph
	vector<pair<int, int>> graph_edges = {{1,2}, {1,3},{1,4},{1,5},{1,6},{2,8},{2,9},{3,7},{4,8},{4,9},{5,9},{6,8},{7,10},{8,10},{9,10}};
	
	
	//SAMPLE EXAMPLE 1
	/*
	int t=8, c=3, edge_no=9; 
	int T_local[t][c] = {{9,7,5}, {8,6,5},{6,5,4},{7,5,3},{5,4,2},{7,6,4},{8,5,3},{6,4,2}};
	vector<pair<int, int>> graph_edges = {{1,4},{2,4},{3,4},{4,5},{4,6},{4,7},{5,8},{6,8},{7,8}};
	*/
	
	// SAMPLE EXAMPLE 2
	/*
	int t=8, c=3, edge_no=10; 
	int T_local[t][c] = {{9,7,5}, {8,6,5},{6,5,4},{7,5,3},{5,4,2},{7,6,4},{8,5,3},{6,4,2},{5,3,2},{7,4,2}};
	vector<pair<int, int>> graph_edges = {{1,2},{1,3},{2,4},{2,5},{3,5},{3,6},{4,7},{5,7},{5,8},{6,8}};
	*/
	
	cout<<"Graph Edges are : "<<endl;
	printVectorElements(graph_edges);
	
	float Pk[3] = {1,2,4}; // Power consumption on core 1, core 2 and core 3 respectively
	float Ps = 0.5; // Power consumption for wireless send
	
	// Creating the graph
	TaskGraph *taskgraph = new TaskGraph(graph_edges, edge_no, t, c);
	taskgraph->displayGraphDetails();
	
	// STEP 1 - INITIAL SCHEDULING
	clock_t begin_initial_scheduling = clock();
	taskgraph->primary_assignment((int *)T_local, T_cloud, t, c);
	taskgraph->task_prioritizing(t-1);
	taskgraph->execution_unit_selection((int *)T_local, T_cloud, Pk, Ps);
	clock_t stop_initial_scheduling = clock();
	double time_initial_scheduling = double(stop_initial_scheduling - begin_initial_scheduling) / CLOCKS_PER_SEC;
	
	// STEP 2 - TASK MIGRATION
	clock_t begin_task_migration = clock();
	taskgraph->Task_Migration((int *)T_local, T_cloud, Pk, Ps);
	clock_t stop_task_migration = clock();
	double time_task_migration = double(stop_task_migration - begin_task_migration) / CLOCKS_PER_SEC;
	
	// Running Time
	cout<<"RUNNING TIME RESULTS : "<<endl;
	cout<<"The Initial Scheduling Phase Runs for Time = "<<time_initial_scheduling<<" seconds."<<endl;
	cout<<"The Task Migration Phase Runs for Time = "<<time_task_migration<<" seconds."<<endl;
	cout<<"The Initial Scheduling Phase Runs for Time = "<<time_initial_scheduling+time_task_migration<<" seconds."<<endl;
	
	return 0;	
}
