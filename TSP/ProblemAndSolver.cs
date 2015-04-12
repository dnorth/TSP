using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using Priority_Queue;

namespace TSP
{
    public class State : PriorityQueueNode
    {
        public double[,] state;
        public double lower_bound;
        public double cost;
        public int cities_left;
        public SortedList edges;
        //public SortedList<int, int> edges;
        //public int[,] edges;
        public State(double[,] state, double lower_bound, double cost, int cities_left, SortedList edges)
        {
            this.state = state;
            this.lower_bound = lower_bound;
            this.cost = cost;
            this.cities_left = cities_left;
            this.edges = edges;
        }
        public State(double[,] state, double lower_bound)
        {
            this.state = state;
            this.lower_bound = lower_bound;
            this.cost = -1;
            this.cities_left = -1;
            this.edges = null;
        }
    }
    class ProblemAndSolver
    {
        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// you are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your node data structure and search algorithm. 
            /// </summary>
            public ArrayList 
                Route;

            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }


            /// <summary>
            ///  compute the cost of the current route.  does not check that the route is complete, btw.
            /// assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here; 
                double cost = 0D;
                
                for (x = 0; x < Route.Count-1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }
                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost; 
            }
        }

        #region private members
        private const int DEFAULT_SIZE = 25;
        
        private const int CITY_ICON_SIZE = 5;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;
        #endregion

        #region public members.
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        public const int DEFAULT_SEED = -1;

        #region Constructors
        public ProblemAndSolver()
        {
            initialize(DEFAULT_SEED, DEFAULT_SIZE);
        }

        public ProblemAndSolver(int seed)
        {
            initialize(seed, DEFAULT_SIZE);
        }

        public ProblemAndSolver(int seed, int size)
        {
            initialize(seed, size);
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// reset the problem instance. 
        /// </summary>
        private void resetData()
        {
            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null; 

            for (int i = 0; i < _size; i++)
                Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.LightGray,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        private void initialize(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            if (seed != DEFAULT_SEED)
                this.rnd = new Random(seed);
            else
                this.rnd = new Random();
            this.resetData();
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size)
        {
            this._size = size;
            resetData(); 
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-15F;
            Font labelFont = new Font("Arial", 10);

            g.DrawString("n(c) means this node is the nth node in the current solution and incurs cost c to travel to the next node.", labelFont, cityBrushStartStyle, new PointF(0F, 0F)); 

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        double BSSF;
        double minCost;
        double currCost;
        int minIndex;
        double[,] initialState;
        HeapPriorityQueue<State> Agenda;
        DateTime end;
        State best_state;

        /// <summary>
        ///  solve the problem.  This is the entry point for the solver when the run button is clicked
        /// right now it just picks a simple solution. 
        /// </summary>
        public void solveProblem()
        {
            DateTime start = DateTime.Now;
            end = start.AddSeconds(60);

            /******* STEP 1 Create the initial State *******/
            init_state();

            /******* STEP 2 Get a 'Quick' BSSF *******/
            quick_solution();
            double init_bound = bound(initialState, 0).lower_bound;
            SortedList edges = new SortedList();
            best_state = new State(initialState, init_bound, 0, Cities.Length, edges);

            /******* STEP 3 Create a new Agenda *******/
            Agenda = new HeapPriorityQueue<State>(1000000);
            Agenda.Clear();

            /******* STEP 4 Add our first Route and BSSF to the agenda *******/
            Agenda.Enqueue(best_state, init_bound);

            /******* STEP 5 Run the algorithm while the Agenda size is greater than 0, We are not over our minute mark, and the BSSF is not equal to the state's lower_bound *******/
            while (Agenda.Count > 0 && DateTime.Now < end && BSSF != Agenda.First.lower_bound)
            {
                //Remove the first state off of the Agenda
                State u = Agenda.Dequeue();

                //Check if the lower bound is less than the current BSSF value
                if (u.lower_bound < BSSF)
                {
                    //As long as the lower bound is less, we will create successors
                    HashSet<State> children = successors(u);
                    if (children == null)
                        continue;
                    //For each successor that we have, check if we ran out of time and check if its lower bound is also less than the BSSF
                    foreach (State child in children)
                    {
                        if (DateTime.Now >= end)
                        {
                            break;
                        }
                        if (child.lower_bound < BSSF)
                        {
                            //As long as the child's lower bound is less than the BSSF, we will check to see if the criterion is met.
                            if (criterion(child))
                            {
                                //If the criterion is met, replace our BSSF and say this is our new best state
                                best_state = child;
                                BSSF = child.cost;
                            }
                            else
                            {
                                //If the criterion is not met, add this state to our agenda but don't just use the lower bound as its priority; we calculate the priority based on a heuristic function
                                Agenda.Enqueue(child, heuristic(child.lower_bound, child.cities_left));
                            }
                        }
                    }             
                }
            }

            //Check to see if our best state is not the initial BSSF we found
            if (best_state.cost != 0)
            {
                Route.Clear();

                //Route.Add(Cities[best_state.edges.Keys[0]]);
                int first_city = (int)best_state.edges.GetKey(0);
                int current_city = first_city;

                //Order the cities for routing and drawing
                while (true)
                {
                    Route.Add(Cities[current_city]);
                    for (int x = 0; x < best_state.edges.Count; x++)
                    {
                        if ((int)best_state.edges.GetKey(x) == current_city)
                        {
                            current_city = (int)best_state.edges.GetByIndex(x);
                            break;
                        }
                    }
                    if (current_city == first_city)
                    {
                        break;
                    }
                }

                bssf = new TSPSolution(Route);
            }

            Program.MainForm.tbCostOfTour.Text = " " + BSSF;

            Program.MainForm.tbElapsedTime.Text = " " + (DateTime.Now - start);

            Program.MainForm.Invalidate();


        }

        // heuristic(lower_bound, cities_left) evaluate heuristic function according to state
        public double heuristic(double lower_bound, int cities_left)
        {
            if (cities_left < 1) return lower_bound;
            return lower_bound + (cities_left * 7);
        }

        // criterion(state) evaluates state to determine if criterion or not
        public bool criterion(State state)
        {
            if (state.edges.Count == Cities.Length) return true;
            return false;
        }

        // evaluate_edge(i,j,state) decide whether to include/exclude, then add to agenda
        public HashSet<State> evaluate_edge(int i, int j, State state)
        {
            
            double[,] include = (double[,])state.state.Clone();
            double[,] exclude = (double[,])state.state.Clone();
            HashSet<State> childrenList = new HashSet<State>();

            //Set the exclude value to inifinity
            exclude[i, j] = double.MaxValue;

            //Get all the values in the columns and rows for the include and set them to infinity
            for (int x = 0; x < Cities.Length; x++)
            {
                include[x, j] = double.MaxValue;
                include[i, x] = double.MaxValue;
            }
            //also set the from-to value to infinity, this one won't be used in an include situation
            include[j, i] = double.MaxValue;

            //If there is more than 1 city, it is possible to exclude that city; if not, then it is only possible to include that city
            if(state.cities_left > 1)
            {
                //Reduce both of their matrices
                State excludeState = bound(exclude, state.lower_bound);
                excludeState.cities_left = state.cities_left;
                excludeState.edges = (SortedList)state.edges.Clone();
                excludeState.cost = state.cost;
                childrenList.Add(excludeState);
            }

            State includeState = bound(include, state.lower_bound);
            includeState.cities_left = state.cities_left - 1;
            includeState.edges = (SortedList)state.edges.Clone();
            includeState.edges.Add(i, j);
            includeState.cost = state.cost + Cities[i].costToGetTo(Cities[j]);
            childrenList.Add(includeState);


            return childrenList;
        }

        // successors(State) find successors for given State
        //For every row and every column in our state matrix, find the first 0 (There will always be a 0 because of reduction) 
        //If we find a 0, make sure we wouldn't be creating a premature cycle by including it
        public HashSet<State> successors(State state)
        {
            for (int x = 0; x < Cities.Length; x++)
            {
                for (int y = 0; y < Cities.Length; y++)
                {
                    if (state.state[x,y] == 0)
                    {
                        if (!premature(x, y, state))
                        {
                            return evaluate_edge(x, y, state);
                        }
                    }
                }
            }
            return null;
        }

        // premature(x, y, state) checks to see if given edge could be used to complete a premature cycle
        // First we make sure there is at least more than 1 city left; if not, we are not in a premature cycle
        // Then we make sure the to city: y never loops around to the from city: x. If that was the case, we would have a premature cycle
        public bool premature(int x, int y, State state)
        {
            if (state.cities_left > 1)
            {
                int first_city = x;
                int current_city = y;
                while(true)
                {
                    if (current_city == first_city)
                    {
                        return true;
                    }
                    if (state.edges.ContainsKey(current_city))
                    {
                        for (int i = 0; i < state.edges.Count; i++)
                        {
                            if ((int)state.edges.GetKey(i) == current_city)
                            {
                                current_city = (int)state.edges.GetByIndex(i);
                                break;
                            }
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            return false;
        }

        // bound(costMatrix) find lower bound for given cost matrix
        //Reduce all of the rows and columns; add the reduction to the lower bound
        public State bound(double[,] state, double bounding)
        {
            double lower_bound = bounding;
            // reduce by row
            for (int x = 0; x < Cities.Length; x++)
            {
                double lowest = double.MaxValue;
                for (int y = 0; y < Cities.Length; y++)
                {
                    if (state[x,y] < lowest)
                    {
                        lowest = state[x,y];
                    }
                }
                if (lowest != 0 && lowest != double.MaxValue)
                {
                    lower_bound += lowest;
                    for (int y = 0; y < Cities.Length; y++)
                    {
                        state[x,y] -= lowest;
                    }
                }
            }

            // reduce by column
            for (int y = 0; y < Cities.Length; y++)
            {
                double lowest = double.MaxValue;
                for (int x = 0; x < Cities.Length; x++)
                {
                    if (state[x,y] < lowest)
                    {
                        lowest = state[x, y];
                    }
                }
                if (lowest != 0 && lowest != double.MaxValue)
                {
                    lower_bound += lowest;
                    for (int x = 0; x < Cities.Length; x++)
                    {
                        state[x, y] -= lowest;
                    }
                }
            }
            return new State(state, lower_bound);
        }

        // quick_solution, not using initial state cost matrix
        public void quick_solution()
        {
            int x;
            int currIndex = 0;
            Route = new ArrayList();

            //Get the next neighbor greedily based on lowest cost
            //While our Route list size is less than the number of cities (meaning we finish once we get to all the cities)
            while (Route.Count < Cities.Length)
            {
                //Get an arbitrarily large minimum cost
                minCost = double.MaxValue;

                //Loop through all the cities to find the one with the least cost
                for (x = 0; x < Cities.Length; x++)
                {
                    //Make sure we're not trying to check against ourself, because that would definitely be the minimum cost and kind of ruin this whole thing.....
                    if (x != currIndex)
                    {
                        //Make sure the Route list doesn't already contain the city we're trying to add. If it does, it won't be a tour so we don't want to use it again
                        if (!Route.Contains(Cities[x]))
                        {
                            //Get the cost from the current index to the next index. If it's less than what we have, make the 'minCost' equal to our 'currCost' and choose our new minimum index
                            currCost = Cities[currIndex].costToGetTo(Cities[x]);
                            if (currCost < minCost)
                            {
                                minCost = currCost;
                                minIndex = x;
                            }
                        }
                    }
                }
                //Once we found the city with the lowest cost that is still not in the tour, we have to add that city to our Route
                currIndex = minIndex;
                Route.Add(Cities[currIndex]);
            }

            //Once we have a complete tour, we need to calculate the cost of it and make it our initial BSSF
            bssf = new TSPSolution(Route);
            BSSF = bssf.costOfRoute();
        }

        // init_state() create an initial matrix based on the costs of all the cities
        // If we are on the diagonal of a matrix, set its value to infinity. Else, find the cost that it takes to go from that row to that column and add it to the matrix
        public void init_state()
        {
            initialState = new double[Cities.Length, Cities.Length];
            for(int i = 0; i < Cities.Length; i++)
            {
                for (int j = 0; j < Cities.Length; j++)
                {
                    if(i == j)
                    {
                        initialState[i, j] = double.MaxValue;
                    }
                    else
                    {
                        initialState[i, j] = Cities[i].costToGetTo(Cities[j]);
                    }
                }
            }    
        }
            
        #endregion
    }
}
