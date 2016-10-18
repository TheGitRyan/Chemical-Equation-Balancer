using System;
using System.Collections;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ChemicalEquationBalancer
{
    
    public class Matrix
    {
        private double[,] entries;

        public double[,] Entries
        {
            get
            {
                return entries;
            }
            set
            {
                entries = value;
            }
        }

        public void SetEntry(int i, int j, double entry)
        {
            entries[i, j] = entry;
        }

        public double GetEntry(int i, int j)
        {
            return entries[i, j];
        }

        public double[] GetRow(int i)
        {
            double[] a = new double[this.Entries.GetLength(1)];

            for (int j = 0; j < this.Entries.GetLength(1); j++)
            {
                a[j] = this.GetEntry(i, j);
            }

            return a;
        }

        public double[] GetColumn(int j)
        {
            double[] a = new double[this.Entries.GetLength(0)];

            for (int i = 0; i < this.Entries.GetLength(0); i++)
            {
                a[i] = this.GetEntry(i, j);
            }

            return a;
        }

        //constructor
        public Matrix(int i, int j)
        {
            this.entries = new double[i, j];
        }

        public Boolean isSquare()
        {
            return entries.GetLength(0) == entries.GetLength(1);
        }

        public Matrix GetInverse()
        {
            double[,] a = entries;
            Matrix inv = this;

            if (a.GetLength(0) != 2) {
                double scal = (1 / this.GetDeterminant());
                inv = this.GetMatrixAdjunct();
                inv = inv.Multiply(scal);
            }
            else
            {
                double scal = (1 / DeterminantOf2X2(this));

                double temp = inv.GetEntry(0, 0);
                inv.SetEntry(0, 0, inv.GetEntry(1,1));
                inv.SetEntry(1,1, temp);
                
                inv.SetEntry(0, 1, (-1) * inv.GetEntry(0,1));
                inv.SetEntry(1, 0, (-1) * inv.GetEntry(1,0));

                inv = inv.Multiply(scal);
            }

            return inv;
        }

        public double GetDeterminant()
        {
            double determinant = 0;
            if (this.isSquare())
            {                
                if (this.Entries.GetLength(0) != 2)
                {
                    for (int i = 0; i < this.Entries.GetLength(0); i++)
                    {
                        Matrix minor = this.GetMinorMatrix(0, i);

                        if (minor.Entries.GetLength(0) == 2)
                        {
                            determinant += (Math.Pow(-1, i) * (DeterminantOf2X2(minor) * this.GetEntry(0, i)));
                        }
                        else
                        {
                            determinant += (Math.Pow(-1, i) * (minor.GetDeterminant() * this.GetEntry(0, i)));
                        }
                    }
                }
                else
                {
                    determinant += DeterminantOf2X2(this);
                }
            }

            return determinant;
        }

        public static double DeterminantOf2X2(Matrix m)
        {
            double[,] a = m.Entries;
            return (a[0, 0] * a[1, 1]) - (a[0, 1] * a[1, 0]);
        }

        public Matrix GetMatrixTranspose()
        {
            double[,] a = entries;
            Matrix tra = new Matrix(a.GetLength(0), a.GetLength(1));
            
            for(int i = 0; i < a.GetLength(0); i++)
            {
                for(int j = 0; j < a.GetLength(1); j++)
                {
                    tra.SetEntry(j, i, this.GetEntry(i, j));
                }
            }

            return tra;
        }

        public Matrix GetMatrixAdjunct()
        {
            Matrix adj = this.GetMatrixOfCofactors();
            return adj.GetMatrixTranspose();
        }

        public Matrix GetMatrixOfCofactors()
        {
            double[,] a = entries;
            Matrix cof = this.GetMatrixOfMinors();

            for(int i = 0; i < a.GetLength(0); i++)
            {
                for(int j = 0; j < a.GetLength(1); j++)
                {
                    if((i + j) % 2 == 1)
                    {
                        cof.SetEntry(i, j, ((-1) * cof.GetEntry(i, j)));
                    }
                }
            }

            return cof;
        }

        //returns a matrix that contains all the minors for a matrix in a 
        public Matrix GetMatrixOfMinors()
        {
            double[,] a = this.Entries;
            Matrix minorMatrix = new Matrix(a.GetLength(0), a.GetLength(1));

            for(int i = 0; i < a.GetLength(0); i++)
            {
                for(int j = 0; j < a.GetLength(1); j++)
                {
                    minorMatrix.SetEntry(i, j, this.GetMinor(i, j));
                }
            }

            return minorMatrix;
        }
            
        //gets a specific minor for the indices that are passed in
        public double GetMinor(int jRemove, int kRemove)
        {
            //creates a new matrix as template for minor, will be added to arraylist that contains all minors
            Matrix m = this.GetMinorMatrix(jRemove, kRemove);
            return m.GetDeterminant();
        }

        //gets a specific minor for the indices that are passed in
        public Matrix GetMinorMatrix(int jRemove, int kRemove)
        {
            //creates a new matrix as template for minor, will be added to arraylist that contains all minors
            Matrix m = new Matrix(this.entries.GetLength(0) - 1, this.entries.GetLength(1) - 1);

            int jMinor = 0;

            //initialize generic j index which will loop through the set of ROWS (particular column) in THIS matrix
            for (int j = 0; j < this.entries.GetLength(0) && jMinor < m.Entries.GetLength(0); j++)
            {

                int kMinor = 0;

                //don't want to loop through that column if it is the one that is removed on this pass
                if (j != jRemove)
                {
                    //if this is not a removed column then an index for looping through the row is initialized
                    for (int k = 0; k < this.entries.GetLength(1) && kMinor < m.Entries.GetLength(1); k++)
                    {

                        //check to make sure this is not a removed row
                        if (k != kRemove)
                        {
                            //if this is not a removed row or column then this enry is valid and should be
                            //entered into minor matrix
                            m.SetEntry(jMinor, kMinor, this.GetEntry(j, k));
                            kMinor++;
                        }
                    }
                    jMinor++;
                }
            }
            if (kRemove == m.entries.Length - 1)
            {
                kRemove = 0;
                jRemove++;
            }
            else {
                kRemove++;
            }

            return m;
        }

        //multiplies two matrices
        public Matrix Multiply(Matrix m)
        {
            double[,] a = this.Entries;
            double[,] b = m.Entries;

            Matrix multi = new Matrix(a.GetLength(0), b.GetLength(1));

            if(a.GetLength(1) == b.GetLength(0))
            {
                for(int i = 0; i < a.GetLength(0); i++)
                {
                    for (int j = 0; j < b.GetLength(1); j++)
                    {
                        double[] row = this.GetRow(i);
                        double[] col = m.GetColumn(j);

                        multi.SetEntry(i, j, LinComb(row, col));
                    }
                }
            }

            return multi;
        }

        //returns linear combination of a andb, they must be the same length
        public static double LinComb(double[] a, double[] b)
        {
            double sum = 0.0;

            for(int i = 0; i < a.Length; i++)
            {
                sum += (a[i] * b[i]);
            }

            return sum;
        }

        //adds rows or columns to make array of entries into a square matrix, fills the added rows with value
        //passed in
        public Matrix ConvertToSquare(double fill)
        {
            int numLinesToAdd = entries.GetLength(0) - entries.GetLength(1);

            //replacement int array that will be loaded in place of current entries;
            double[,] squareReplace;
            if (numLinesToAdd < 0) {
                //if numLinesToAdd is negative then there are too many columns (not enough rows), so lines are added
                //to rows
                numLinesToAdd *= -1;
                squareReplace = new double[entries.GetLength(0) + numLinesToAdd, entries.GetLength(1)];
            }
            else if(numLinesToAdd > 0)
            {
                //if numLinesToAdd is positive then there are too many rows (not enough columns), so lines are added 
                //to columns
                squareReplace = new double[entries.GetLength(0), entries.GetLength(1) + numLinesToAdd];
            } else
            {
                return this;
            }

            //fills the replacement array with value from entries, or if outside of the range of entries, then fills
            //with fill
            for(int i = 0; i < squareReplace.GetLength(0); i++)
            {
                for (int j = 0; j < squareReplace.GetLength(1); j++)
                {
                    if (i < entries.GetLength(0) && j < entries.GetLength(1)) {
                        squareReplace[i, j] = entries[i, j];
                    } 
                    else
                    {
                        squareReplace[i, j] = fill;
                    }
                }

            }

            entries = squareReplace;
            return this;
        }

        //overloaded ConvertToSquare with no arguement, autoloads with 0
        public Matrix ConvertToSquare()
        {
            return this.ConvertToSquare(0);
        }

        //multiply overloaded for scalar
        public Matrix Multiply(double n)
        {
            double[,] a = this.Entries;
            Matrix m = new Matrix(a.GetLength(0), a.GetLength(1));

            for(int i = 0; i < a.GetLength(0); i++)
            {
                for(int j = 0; j < a.GetLength(1); j++)
                {
                    m.SetEntry(i, j, (n * this.GetEntry(i, j)));
                }
            }

            return m;
        }

        public override string ToString()
        {
            string s = "[";

            for(int i = 0; i < this.Entries.GetLength(0); i++)
            {
                if(i == 0)
                {
                    s += "[";
                }
                else
                {
                    s += " [";
                }
                
                for(int j = 0; j < this.Entries.GetLength(1); j++)
                {
                    s += this.GetEntry(i , j);
                    if(j != this.Entries.GetLength(1) - 1)
                    {
                        s += ", ";
                    }
                }
                s += "]";

                if(i != this.Entries.GetLength(0) - 1)
                {
                    s += "\n";
                }
            }
            s += "]";

            return s;
        }
    }
}
