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
        public int[,] route;
        public State(double[,] state, double lower_bound, double cost, int cities_left, int[,] route)
        {
            this.state = state;
            this.lower_bound = lower_bound;
            this.cost = cost;
            this.cities_left = cities_left;
            this.route = route;
        }
        public State(double[,] state, double lower_bound)
        {
            this.state = state;
            this.lower_bound = lower_bound;
            this.cost = -1;
            this.cities_left = -1;
            this.route = null;
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
            int[,] route = new int[Cities.Length, 2];
            best_state = new State(initialState, init_bound, BSSF, Cities.Length - 1, route);

            /******* STEP 3 Create a new Agenda *******/
            // TODO: Agenda will be our Priority Queue (do we need to implement the multiplication thing when we add?? I dunno yet.)
            Agenda = new HeapPriorityQueue<State>(1000);
            Agenda.Clear();
            Agenda.Enqueue(best_state, init_bound);
            while (Agenda.Count > 0 && DateTime.Now < end && BSSF != Agenda.First.lower_bound)
            {
                State u = Agenda.Dequeue();
                if (u.lower_bound < BSSF)
                {
                    HashSet<State> children = successors(u);
                    if (children == null)
                    {
                        String str = "OH NO";
                    }
                    foreach(State child in children)
                    {

                        child.cost = heuristic(child.lower_bound, child.cities_left);
                        if (DateTime.Now >= end)
                        {
                            break;
                        }
                        //This is the weird requirement we're not sure of yet. :D
                        if (child.lower_bound < BSSF)
                        {
                            if(criterion(child.state))
                            {
                                BSSF = child.cost;
                            }
                            else
                            {
                                Agenda.Enqueue(child, child.cost);
                            }
                        }
                    }
                        
                }
            }

            Route.Clear();

            int[] ordered_route = get_ordered_route();

            for (int i = 0; i < ordered_route.Length; i++)
            {
                Route.Add(Cities[ordered_route[i]]);
            }

            bssf = new TSPSolution(Route);

            Program.MainForm.tbCostOfTour.Text = " " + BSSF;

            Program.MainForm.tbElapsedTime.Text = " " + (DateTime.Now - start);

            Program.MainForm.Invalidate();
            /******* STEP 4 Add our first Route and BSSF to the agenda *******/
            // TODO: We will add our state and the BSSF to the agenda

            /******* STEP 5 Add our first Route and BSSF to the agenda *******/

        }

        // get_ordered_route() returns ordered route of cities after BSSF is determined
        public int[] get_ordered_route()
        {
            int[] ordered_route = new int[Cities.Length];
            int[,] state_route = best_state.route;

            ordered_route[0] = state_route[0, 0];
            ordered_route[1] = state_route[0, 1];
            state_route[0, 0] = int.MaxValue;
            state_route[0, 1] = int.MaxValue;

            int current_city = ordered_route[1];

            for (int i = 2; i < Cities.Length; i++)
            {
                int next_city = -1;
                for (int x = 0; x < Cities.Length; x++)
                {
                    if (state_route[x, 0] == current_city) next_city = state_route[x, 1];
                    else if (state_route[x, 1] == current_city) next_city = state_route[x, 0];
                    if (next_city != -1)
                    {
                        state_route[x, 0] = int.MaxValue;
                        state_route[x, 1] = int.MaxValue;
                        break;
                    }
                }
                ordered_route[i] = next_city;
                current_city = next_city;
            }
            return ordered_route;
        }

        // heuristic(lower_bound, cities_left) evaluate heuristic function according to state
        public double heuristic(double lower_bound, int cities_left)
        {
            if (cities_left < 1) return lower_bound;
            return lower_bound + (cities_left * 7);
        }

        // criterion(state) evaluates state to determine if criterion or not
        public bool criterion(double[,] state)
        {
            for (int x = 0; x < Cities.Length; x++)
            {
                for (int y = 0; y < Cities.Length; y++)
                {
                    if (state[x,y] != double.MaxValue)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        // evaluate_edge(i,j,state) decide whether to include/exclude, then add to agenda
        public HashSet<State> evaluate_edge(int i, int j, State state)
        {
            double[,] include = state.state;
            double[,] exclude = state.state;
            HashSet<State> childrenList = new HashSet<State>();

            //Set the exclude value to inifinity
            exclude[i, j] = double.MaxValue;

            //Get all the values in the columns and rows for the include and set them to infinity
            for (int x = 0; x < Cities.Length; x++)
            {
                include[x, j] = double.MaxValue;
                include[i, x] = double.MaxValue;
            }
            include[j, i] = double.MaxValue;
            //Reduce both of their matrices
            State includeState = bound(include, state.lower_bound);
            State excludeState = bound(exclude, state.lower_bound);

            includeState.cities_left = state.cities_left - 1;
            excludeState.cities_left = state.cities_left - 1;
            includeState.route = state.route;
            excludeState.route = state.route;
            includeState.route[Cities.Length - 1 - state.cities_left, 0] = i;
            includeState.route[Cities.Length - 1 - state.cities_left, 1] = j;
            includeState.cost = heuristic(includeState.lower_bound, includeState.cities_left);
            excludeState.cost = heuristic(excludeState.lower_bound, excludeState.cities_left);

            childrenList.Add(includeState);
            childrenList.Add(excludeState);


            return childrenList;
        }
            /*
            double exclude_bound = exclude_state.lower_bound;
            double include_bound = include_state.lower_bound;
            double lowest_bound = 0;
            if (exclude_bound < include_bound)
            {
                lowest_bound = exclude_bound;
                lowest_state = exclude_state.state;
            }
            else
            {
                lowest_bound = include_bound;
                lowest_state = include_state.state;
            }*/
            /*
        }*/

        // successors(State) find successors for given State
        public HashSet<State> successors(State state)
        {
            for (int x = 0; x < Cities.Length; x++)
            {
                for (int y = 0; y < Cities.Length; y++)
                {
                    if (state.state[x,y] != double.MaxValue)
                    {
                        return evaluate_edge(x, y, state);
                    }
                }
            }
            return null;
        }

        // bound(costMatrix) find lower bound for given cost matrix
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

        // return factorial of given number
        public int factorial(int number)
        {
            int fact = number;
            for (int x = 1; x < number; x++)
            {
                fact *= x;
            }
            return fact;
        }

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
