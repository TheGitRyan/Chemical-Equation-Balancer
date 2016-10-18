using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace ChemicalEquationBalancer
{
    class Program
    {
        static void Main(string[] args)
        {
            Boolean run = true;
            while (run)
            {
                Console.Write("Enter an equation: ");
                string eq = Console.ReadLine();
                if (eq.ToLower().Equals("quit"))
                {
                    run = false;
                }
                else
                {
                    //get a list of all elements and chemicals
                    ArrayList elements = AllUniqueElements(eq);
                    ArrayList chemicals = AllUniqueChemicals(eq);

                    //get matrix with for all chems but last and matrix for last chem
                    Matrix A = PopulateChemMatrix((string[])chemicals.ToArray(typeof(string)), (string[])elements.ToArray(typeof(string)), eq);
                    Matrix B = PopulateEquivalentMatrix((string[])chemicals.ToArray(typeof(string)), (string[])elements.ToArray(typeof(string)));

                    //calculate coefficients
                    Matrix C = (A.GetInverse().Multiply(B)).Multiply(A.GetDeterminant());

                    //make and load an arrayList with all values in Matrix c
                    double[] coefficients = new double[C.GetColumn(0).Length];
                    Array.Copy(C.GetColumn(0), coefficients, C.GetColumn(0).Length);

                    double last = coefficients[coefficients.Length - 1];
                    if (coefficients.Length == ((string[])chemicals.ToArray(typeof(string))).Length)
                    {
                        if (last != A.GetDeterminant())
                        {
                            coefficients[coefficients.Length - 1] = A.GetDeterminant();
                        }
                    } else
                    {
                        //if the array of coefficients is one short then add the 
                        //det of A to the end as the last coefficient

                        ArrayList c = new ArrayList(C.GetColumn(0));
                        c.Add(A.GetDeterminant());
                        coefficients = (double[])c.ToArray(typeof(double));
                    }

                    coefficients = MakeAllPositive(coefficients); //all coefficients must be positive

                    coefficients = ReduceSet(coefficients); //reduces all numebrs into lowest terms

                    string final = FillCoefficients(coefficients, eq);
                    
                    Console.WriteLine(final);
                }
            }
        }

        //reduce set to lowest terms
        public static double[] ReduceSet(double[] dArr)
        {
            for(int i = 0; i < dArr.Length; i++)
            {
                double currentD = dArr[i];
                for(int j = 0; j < dArr.Length; j++)
                {
                    double checkD = dArr[j];
                    double greatestCommon = GCD(currentD, checkD);
                    if(DividesAll(greatestCommon, dArr))
                    {
                        DivideAll(greatestCommon, dArr);
                    }
                }
            }
            return dArr;
        }

        public static Boolean DividesAll(double div, double[] dArr)
        {
            foreach (double d in dArr)
            {
                if (d % div != 0)
                {
                    return false;
                }
            }
                return true;
        } 

        public static double[] DivideAll(double div, double[] dArr)
        {
            for(int i = 0; i < dArr.Length; i++)
            {
                dArr[i] = (dArr[i] / div);
            }
            return dArr;
        }

        public static double GCD(double a, double b)
        {
            while (a != 0 && b != 0)
            {
                if (a > b)
                    a %= b;
                else
                    b %= a;
            }

            if (a == 0)
                return b;
            else
                return a;
        }
        
        //returns a string with the coefficients of theequation added to balance the equation
        public static string FillCoefficients(double[] coef, string eq)
        {
            string[] eqParse = eq.Split(' ');

            string balance = "";

            int indexCoef = 0;
            int indexEq = 0;
            while (indexCoef < coef.Length)
            {
                int caseNum = 0;
                if (eqParse[indexEq].Equals("->") || eqParse[indexEq].Equals("+"))
                {
                    caseNum = 1;
                } else
                {
                    caseNum = 2;
                }

                switch (caseNum)
                {
                    case 1:
                        balance += " " + eqParse[indexEq] + " ";
                        indexEq++;
                        break;
                    case 2:
                        balance += "(" + coef[indexCoef] + ")" + eqParse[indexEq];
                        indexCoef++;
                        indexEq++;
                        break;
                    default:
                        break;
                }
            }

            return balance;
        }


        //returns whether or not a chemical is a reactant 
        public static Boolean IsReactant(string chem, string eq)
        {
            string[] eqSplit = eq.Split(' ');

            //get index of the yield sign
            int indexYield = -1;
            for (int i = 0; i < eqSplit.Length; i++)
            {
                if (eqSplit[i].Equals("->"))
                {
                    indexYield = i;
                    break;
                }
            }

            int indexChem = Array.BinarySearch(eqSplit, chem);

            return indexChem < indexYield;
        }

        //returns whether or not a chemical is a product
        public static Boolean IsProduct(string chem, string eq)
        {
            string[] eqSplit = eq.Split(' ');

            //get index of the yield sign
            int indexYield = -1;
            int indexChem = -1;
            for (int i = 0; i < eqSplit.Length; i++)
            {
                if (eqSplit[i].Equals("->"))
                {
                    indexYield = i;
                    break;
                }
            }

            for (int i = 0; i < eqSplit.Length; i++)
            {
                if (eqSplit[i].Equals(chem))
                {
                    indexChem = i;
                    break;
                }
            }

            return indexChem > indexYield;
        }

        public static ArrayList AllUniqueChemicals(string unsplit)
        {
            string[] unparsed = unsplit.Split(' ');

            ArrayList allUniqueChemicals = new ArrayList();
            foreach (string s in unparsed)
            {
                if (!(allUniqueChemicals.Contains(s) || s.Equals("+") || s.Equals("->")))
                {
                    allUniqueChemicals.Add(s);
                }
            }

            return allUniqueChemicals;
        }

        //takes split but unparsed equation and returns an array of all the unique elements
        public static ArrayList AllUniqueElements(string unsplit)
        {
            string[] unparsed = unsplit.Split(' ');

            ArrayList allUniqueElements = new ArrayList();
            //this loop will loop through each token in eq
            for (int i = 0; i < unparsed.Length; i++)
            {
                string chemical = unparsed[i];

                if (chemical != "+" && chemical != "->")
                {
                    //this loop will go through each token (chemical or operator) and break it up into elements
                    //check to make sure there is a character after the current one with 
                    //which to do checks in switch statement

                    for (int j = 0; j < chemical.Length; j++)
                    {
                        //boundary case where all rest of chars are numbers, there will be no
                        //more elements, but number will be added to list of elements
                        if (IsNumeric(chemical[j]) && RestAreNumeric(j, chemical))
                        {
                            break;
                        }

                        //skip over numbers and only execute creation of an element
                        //if we find a new Uppercase number
                        while (IsNumeric(chemical[j]) && j < chemical.Length - 1)
                        {
                            j++;
                        }

                        string element = chemical[j] + "";

                        if (j < chemical.Length - 1)
                        {
                            int nextLetterSwitch;
                            char nextChar = chemical[j + 1];

                            if (IsNumeric(nextChar))
                            {
                                nextLetterSwitch = 1;
                            }
                            else if (IsUpperCase(nextChar))
                            {
                                nextLetterSwitch = 2;
                            }
                            else
                            {
                                nextLetterSwitch = 3;
                            }

                            switch (nextLetterSwitch)
                            {
                                case 1: //case where next char is a number
                                    if (!ListContains(allUniqueElements, element))
                                    {
                                        allUniqueElements.Add(element);

                                    }

                                    break;
                                case 2: //case where next char is uppercase
                                    goto case 1;
                                case 3: //case where next char is lowercase
                                    element += nextChar;
                                    j++;
                                    goto case 1;
                                default:
                                    break;

                            }
                        }
                    }
                }
            }

            return allUniqueElements;
        }

        //populates the chemical matrix with the amount of each element contained in each chemical
        public static Matrix PopulateChemMatrix(string[] chems, string[] elmts, string eq) {
            Matrix A = new Matrix(elmts.Length, chems.Length - 1);

            for (int i = 0; i < A.Entries.GetLength(0); i++)
            {
                string element = elmts[i];

                for (int j = 0; j < A.Entries.GetLength(1); j++)
                {
                    string chemical = chems[j];

                    if (chemical.Contains(element))
                    {
                        //gets the character after the element which indicates ther number of atoms,
                        //if it is the beginning of the next element then that element only appears once in
                        //the chemical

                        if (chemical.Equals(element))
                        {
                            A.SetEntry(i, j, 1);
                        }
                        else
                        {
                            int endOfAtom = FindLastIndex(chemical, element) + 1;
                            char num = '\0';
                            //if the atom is only one and it is the last element in the chemical formula
                            //then it will run out of bounds on index, 
                            if (endOfAtom != chemical.Length)
                            {
                                num = chemical[endOfAtom];

                                if (IsUpperCase(num))
                                {
                                    //the character after is an uppercase letter and therefore the previous element appears
                                    //only once

                                    if (IsProduct(chemical, eq))
                                    {
                                        A.SetEntry(i, j, -1);
                                    }
                                    else
                                    {
                                        A.SetEntry(i, j, 1);
                                    }
                                }
                                else
                                {
                                    //character is a number and therefore should be converted from to an int
                                    //ascii character for 0 - 9 is 48 - 57 so must be adjusted
                                    int numInt = ToInt(num);

                                    endOfAtom++;
                                    //while there is still more of the number ot get after the atom it will
                                    //continue to append the digits 
                                    while (endOfAtom < chemical.Length && IsNumeric(chemical[endOfAtom]))
                                    {
                                        numInt = AppendNum(numInt, chemical[endOfAtom]);
                                        endOfAtom++;
                                    }

                                    if (IsProduct(chemical, eq))
                                    {
                                        numInt *= -1;
                                    }

                                    A.SetEntry(i, j, numInt);
                                }
                            }
                            else
                            {
                                if (IsProduct(chemical, eq))
                                {
                                    A.SetEntry(i, j, -1);
                                }
                                else
                                {
                                    A.SetEntry(i, j, 1);
                                }
                            }
                        }
                    }
                    else
                    {
                        A.SetEntry(i, j, 0);
                    }
                }
            }

            A.ConvertToSquare(1);

            return A;
        }

        public static Matrix PopulateEquivalentMatrix(string[] chemicals, string[] elmts)
        {
            Matrix A = new Matrix(elmts.Length, 1);

            string chemical = chemicals[chemicals.Length - 1];

            for (int i = 0; i < elmts.Length; i++)
            {
                string element = elmts[i];


                if (chemical.Contains(element))
                {
                    //gets the character after the element which indicates ther number of atoms,
                    //if it is the beginning of the next element then that element only appears once in
                    //the chemical

                    if (chemical.Equals(element))
                    {
                        A.SetEntry(i, 0, 1);
                    }
                    else
                    {
                        int endOfAtom = FindLastIndex(chemical, element) + 1;
                        char num = '\0';
                        //if the atom is only one and it is the last element in the chemical formula
                        //then it will run out of bounds on index, 
                        if (endOfAtom != chemical.Length)
                        {
                            num = chemical[endOfAtom];

                            if (IsUpperCase(num))
                            {
                                //the character after is an uppercase letter and therefore the previous element appears
                                //only once


                                A.SetEntry(i, 0, 1);

                            }
                            else
                            {
                                //character is a number and therefore should be converted from to an int
                                //ascii character for 0 - 9 is 48 - 57 so must be adjusted
                                int numInt = ToInt(num);

                                endOfAtom++;
                                //while there is still more of the number ot get after the atom it will
                                //continue to append the digits 
                                while (endOfAtom < chemical.Length && IsNumeric(chemical[endOfAtom]))
                                {
                                    numInt = AppendNum(numInt, chemical[endOfAtom]);
                                    endOfAtom++;
                                }


                                A.SetEntry(i, 0, numInt);
                            }
                        }
                        else
                        {

                            A.SetEntry(i, 0, 1);

                        }
                    }
                }
                else
                {
                    A.SetEntry(i, 0, 0);
                }
            }

            return A;
        }

        public static double[] MakeAllPositive(double[] dArr)
        {
            for(int i = 0; i < dArr.Length; i++)
            {
                dArr[i] = Math.Abs(dArr[i]);
            }
            return dArr;
        }

        //searches a string for a substring, and returns the last index of the substring
        public static int FindLastIndex(string str, string subStr)
        {
            //if they are equals then return one less than the length
            if (str.Equals(subStr))
            {
                return str.Length - 1;
            }

            //specific case where the substr is only one character then can loop through every character in chemical
            if(subStr.Length == 1)
            {
                char c = subStr[0];
                for(int i = 0; i < str.Length; i++)
                {
                    if(c == str[i])
                    {
                        return i;
                    }
                }
            }
                        
            for (int i = 0; i < str.Length - subStr.Length + 1; i++)
            {
                string checkString = str.Substring(i, subStr.Length);
                    
                if (checkString.Equals(subStr))
                {
                    return i + subStr.Length - 1;
                }
            }
            return -1;          
        }


        //bunch of small methods to handle cases in above methods

        public static Boolean ListContains(ArrayList AL, String s)
        { 
            foreach(object o in AL)
            {
                string ss = (string)o;
                
                if ((ss.Equals(s)))
                {
                    return true;
                }
            }
            return false;
        }

            //returns whether all characters in the string past the index given 
            //are numbers
        public static Boolean RestAreNumeric(int j, string s)
        {
            for(int i = j; i < s.Length; i++)
            {
                if (!IsNumeric(s[i]))
                {
                    return false;
                }
            }
            return true;
        }

        public static Boolean IsUpperCase(char c)
        {
            return ((c >= 65) && (c <= 90));
        }

        public static Boolean IsNumeric(char c)
        {
            return ((c >= 48) && (c <= 57));
        }

        public static int ToInt(char c)
        {
            return c - 48;
        }

        public static int AppendNum(int i, char c)
        {
            if (IsNumeric(c))
            {
                return (i * 10) + ToInt(c);
            } else {
                return -1;
            }
        }
        
    }

}
