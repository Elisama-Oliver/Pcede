
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iostream>

using namespace std;
#define EPSILON 10e-6
#define RC_EPS 1.0e-6
#define tab '\t'


ILOSTLBEGIN


typedef IloArray<IloExtractableArray> Extractable2;
typedef IloArray<Extractable2> Extractable3;

typedef IloArray<IloRangeArray> IloRangeMatrix2;
typedef IloArray<IloRangeMatrix2> IloRangeMatrix3;

typedef IloArray<IloNumVarArray> NumVarMatrix2;
typedef IloArray<NumVarMatrix2> NumVarMatrix3;

typedef IloArray<IloNumArray> IloNumMatrix2;
typedef IloArray<IloNumMatrix2> NumMatrix3;

static void report1 (IloCplex& cutSolver,
		IloNum& periodos,
		IloArray<IloNumVarArray> x,
        IloArray<IloRangeArray> Fill);

static void report2 (IloAlgorithm& patSolver,
		IloNum& periodos,
		IloArray<IloNumVarArray> Use,
		IloObjective obj);

static void report3 (IloCplex& cutSolver,
		IloNum& periodos,
		IloNum& numItens,
        IloArray<IloNumMatrix2> S,
		IloNumVarArray T,
		IloNumVarArray E,
		IloNumVarArray TP,
		IloNumVarArray TI,
        IloArray<IloNumVarArray>EstoqueItens,
		IloArray<IloNumVarArray> x,
		IloArray<IloNumVarArray> v,
		IloArray<NumVarMatrix2> Y);


static void Arredondamento(IloCplex& cplex,
                           IloNum& periodos,
                           IloNumArray& AtrasoRel,
                           IloNumArray& AdiantamentoRel,
                           IloNumArray& TempoPreparoRel,
                           IloArray<IloNumVarArray>x,
                           IloArray<IloNumArray> sol,
                           IloArray<IloNumArray> setupRel,
                           IloNum objFO);




/// VALOR MAXIMO DE ENTRE DOIS VALORES
int max(int a, int b)
{
	int x = a > b ? a:b;
	return (x);
}
int max(int a, int b);

typedef IloArray<IloNumArray> IloMatrix2;
typedef IloArray<IloNumMatrix2> IloMatrix3;




///========================================================================================================///
///                                         PROGRAMA PRINCIPAL                                             ///
///========================================================================================================///


/* Função principal */
void  principal(char instancia[]){

	IloEnv env;
	try{

		IloModel cutOpt (env);


		IloTimer crono(env);  /// Variável para coletar o tempo
		crono.start();

		IloTimer crono2(env);  /// Variável para coletar o tempo
		crono2.start();

		IloTimer crono3(env);  /// Variável para coletar o tempo
		crono3.start();

		IloTimer cronoNo(env);  /// Variável para coletar o tempo
		cronoNo.start();
		/**********************************************************************************************************************************************/
		/****************************                            LEITURA DE DADOS                          ********************************************/
		/**********************************************************************************************************************************************/

		IloInt i, j, t;
		IloNum CompObjeto, periodos, numItens, CustoSetup;

		ifstream entrada(instancia);
		entrada  >> CompObjeto;
		 cout << " Comprimento do objeto = " << CompObjeto << endl;

		entrada  >> numItens;
		 cout << " Numero de itens = " << numItens << endl;

		entrada  >> periodos;
		 cout << " Numero de periodos = " <<periodos << endl;


		IloNumArray size(env,numItens+1);
		for(int i=1; i<=numItens; i++)
			entrada >> size[i];
		 cout << " comprimento dos itens = " << size << endl;

		IloArray<IloNumArray> demanda(env,periodos+1);
		for(int t=1; t<=periodos; t++){
			demanda[t] = IloNumArray(env, numItens+1);
			for(int i=1; i<=numItens; i++){

					entrada >> demanda[t][i];

			}
		}
		         cout << " demanda = " << demanda << endl;



		IloNumArray data(env,periodos+1);
		for(int t=1; t<=periodos; t++)
			entrada >> data[t];
		 cout << " data = " << data << endl;

		IloNumArray alfa(env,periodos+1);
		for(int t=1; t<=periodos; t++)
			entrada >> alfa[t];
		cout << " alfa = " << alfa << endl;

		IloNumArray beta(env,periodos+1);
		for(int t=1; t<=periodos; t++)
			entrada >> beta[t];
		 cout << " beta = " << beta << endl;

		// beta[t] = 1 + rand() % 10;
		// alfa[t] = 1 + rand() % 10;

	IloArray<IloNumArray> ce(env,periodos+1);
		for(int t=1; t<=periodos; t++){
			ce[t] = IloNumArray(env, numItens+1);
			for(int i=1; i<=numItens; i++){

					entrada >> ce[t][i];

			}
		}
    cout << " CUSTO ESTOQUE = " << ce << endl;


		 IloNum numMaxPadroes = numItens;


    /******************* CORTE DE ESTOQUE ********************/
		IloArray<IloNumVarArray> x(env, periodos+1);
		for(IloInt t=1; t<=periodos; t++)
			x[t] = IloNumVarArray(env);



    /******************* CORTE DE ESTOQUE ********************/
		IloArray<IloNumVarArray> v(env, periodos+1);
		for(IloInt t=1; t<=periodos; t++)
			v[t] = IloNumVarArray(env);

    /******************* CORTE DE ESTOQUE ********************/
		IloArray<NumVarMatrix2> Y(env, periodos + 2); /// Aqui!!!
		for(int t =1; t <=periodos; t++){
			Y[t] = NumVarMatrix2(env,  numMaxPadroes +1);
			for(int j = 0; j <=  numMaxPadroes ; j++){
				Y[t][j] = IloNumVarArray(env, numMaxPadroes +1, 0.0,  1.0, ILOFLOAT);
			}
		}
		Y[periodos + 1] = NumVarMatrix2(env,  1); /// Aqui!!!
		Y[periodos + 1][0] = IloNumVarArray(env, numMaxPadroes +1, 0.0,  1.0, ILOFLOAT);/// Aqui!!!

    /*********************  Custo de tempo padrao  *****************************/
		IloArray<IloNumVarArray> TC(env, periodos+1);
		for(IloInt t=1; t<=periodos; t++){
			TC[t] = IloNumVarArray(env,  numMaxPadroes +1, 0.0,  IloInfinity,ILOFLOAT);
		}


    /*********************  ESTOQUE  *****************************/
		IloArray<IloNumVarArray> EstoqueItens(env, periodos+1);
		for(IloInt t=1; t<=periodos; t++){
			EstoqueItens[t] = IloNumVarArray(env,  numItens+1, 0.0,  IloInfinity, ILOFLOAT);
		}

    /********************* VALORES DE TEMPO *****************************/
		IloNumVarArray T( env,periodos+1, 0.0, IloInfinity, ILOFLOAT);
		IloNumVarArray E( env,periodos+1, 0.0, IloInfinity, ILOFLOAT);
		IloNumVarArray TP(env,periodos+1, 0.0, IloInfinity, ILOFLOAT);
		IloNumVarArray TI(env,periodos+1, 0.0, IloInfinity, ILOFLOAT);
		IloNumVarArray SI(env,periodos+1, 0.0, IloInfinity, ILOFLOAT);



		/*** VALOR B-GRANDE ***/
		IloNum M = 1000000;




          IloArray<IloNumArray> tpcx(env,periodos+1);
		for(int t=1; t<=periodos; t++){
			tpcx[t] = IloNumArray(env, x[t].getSize() +1);
			for(int j=1; j<= x[t].getSize() +1; j++){
				 tpcx[t][j] = 1;

			}
		}
		     //   cout << " beta = " << tpcx << endl;


		  IloArray<IloNumArray> c(env,periodos+1);
		for(int t=1; t<=periodos; t++){
			c[t] = IloNumArray(env, x[t].getSize() +1);
			for(int j=1; j<= x[t].getSize()  ; j++){
				 c[t][j] = 1;

			}
		}



/********************************************************************/
///   ADICIONA COLUNAS INICIAIS : COLUNAS HOMOGENEAS               ///
/********************************************************************/


      IloArray<IloNumArray> PadraoCorte1(env, numItens+1);
        for(IloInt j=0; j<=numItens; j++)
            PadraoCorte1[j] = IloNumArray(env, numItens+1);

        for(IloInt j=0; j<=numItens; j++){
            for(IloInt i=1; i<=numItens; i++){
                if(i==j)
             //        PadraoCorte1[j][i] = 1;
                    PadraoCorte1[j][i] = (int(CompObjeto / size[i]));
                else
                    PadraoCorte1[j][i] = 0;
             }
        }


    IloMatrix3 S(env, periodos+1);
    for (int t = 1; t <= periodos; t++) {
      S[t] = IloMatrix2(env, x[t].getSize() +1);
      for (int j = 0; j <= x[t].getSize() ; j++) {
        S[t][j] = IloNumArray(env, x[t].getSize() +1 );
         for (int k = 1; k <= x[t].getSize() ; k++) {

               int soma = 0;

           //   for (int i = 1; i <= numItens; i++){
           //     soma += abs(PadraoCorte1[k][i] - PadraoCorte1[j][i]);

           //   }

             //  S[t][j][k]  =  soma;
               S[t][j][k]  = 1;

         }
      }
    }

   // for (int t = 1; t <= periodos; t++)
    //     cout << " S = " << S[t] << endl;



/*****************************************************************************************************************************************/

		IloExpr foobj(env);


	//	for (int k = 1; k <=  numMaxPadroes ; k++){
	//		obj += S[1][0][k]*Y[1][0][k];
	//	}


		for (int t = 1; t <= periodos; t++) {

			//for (int j = 1; j <= numMaxPadroes; j++){
			//	for (int k = 1; k <= numMaxPadroes; k++ && j!=k)
			//		obj += S[t][j][k]*Y[t][j][k];
			//}

      for (int i = 1; i <=  numItens; i++)
             foobj += ce[t][i]*EstoqueItens[t][i];

			foobj += alfa[t] * T[t];
			foobj += beta[t] * E[t];

		}



        IloObjective RollsUsed = IloAdd(cutOpt, IloMinimize(env, foobj));






/**********************************************************************************************************************************************/
/****************************                  RESTRIÇÕES DO MODELO COM ILORANGE                        ***************************************/
/**********************************************************************************************************************************************/


     IloArray<IloRangeArray>Fill(env,periodos+1);
         for (IloInt t = 1; t <= periodos ; t++){
          Fill[t] = IloRangeArray(env, numItens+1);
             for(IloInt i=1; i<=numItens; i++){
              if(t==1)
               Fill[t][i] = IloRange(env, 0.0, - demanda[t][i]  - EstoqueItens[t][i], 0.0);
              else
                Fill[t][i] = IloRange(env, 0.0, EstoqueItens[t-1][i] - demanda[t][i] - EstoqueItens[t][i], 0.0);
         cutOpt.add(Fill[t][i]);
           }
        }





    IloRangeMatrix2 FillSetup(env, periodos+1);
     for (IloInt t = 1; t<=periodos; t++){
      FillSetup[t] = IloRangeArray(env,numMaxPadroes+1);
       for(IloInt j=1; j<=numMaxPadroes; j++){
         FillSetup[t][j] = IloRange(env, - IloInfinity, 0.0 );
   cutOpt.add(FillSetup[t][j]);
       }
     }






			//  Restrições de estoque Atraso ou Adiantamento
        for (IloInt t = 0; t < periodos ; t++){
	       cutOpt.add(TI[t] + TP[t] + E[t] - T[t] == data[t]);
        }


			//  Restrições de estoque
        for (IloInt t = 1; t <= periodos ; t++){
            if(t==1)
				cutOpt.add(TI[t] >= 0);
			else
				cutOpt.add(TI[t] >= data[t-1] - E[t-1] + T[t-1]);

        }


			/******               RESTRIÇÃO DE ORDENAMENTO DE PADROES NA SEUQNECIA DE CORTE            *********/
  /*  for (IloInt t = 1; t <= periodos ; t++){
		for (int k = 0; k <= numMaxPadroes; k++) { /// Aqui!!!
				IloExpr soma1(env);
				for (int j = 0; j <=numMaxPadroes; j++) {
					if(j!=k)
						soma1 += Y[t][j][k];
				}
				cutOpt.add(soma1 == v[t][k]);

			}
    }



		for (int t = 1; t <= periodos; t++) { /// Substituição de periodos - 1 por periodos
			for (int j = 1; j <= numMaxPadroes; j++) {
				IloExpr soma3(env);
				for (int k = 1; k <= numMaxPadroes; k++) {
					if(j!=k)
						soma3 += Y[t][j][k];
				}
				cutOpt.add(soma3 + Y[t+1][0][j] == v[t][j]);
			}
		}




		for (int t = 1; t <= periodos-1; t++) {
			for (int j = 1; j <= numMaxPadroes; j++) {
				IloExpr soma4(env);
				IloExpr soma5(env);
				for (int k = 1; k <= numMaxPadroes; k++) {
					if(k!=j){
						soma4 += Y[t][k][j];
						soma5 += Y[t][j][k];
					}
				}
				cutOpt.add( Y[t][0][j] + soma4 == soma5 + Y[t+1][0][j]);

			}
		}





		/*******************************************************************************/
		///   IMPOE DE NO INICIO DE CADA PERIODO DE TEMPO A MAQUINA SEJA CONFIGURADA  ///
		/*******************************************************************************/
	/*	for (int t = 1; t <= periodos; t++) {
			IloExpr soma6(env);
			for (int k = 1; k <= numMaxPadroes; k++) {
				soma6 += Y[t][0][k];
			}
			cutOpt.add( soma6 == 1);
		}



		/*****************************************************************************************/
		///              CONTA O TEMPO DE CORTE DO PRIMEIRO PADRAO DE CORTE DO PERIODO  T=1     /////
		/*****************************************************************************************/
 /*       for (int t = 1; t <= periodos; t++) {
           for (int k = 1; k <= numMaxPadroes; k++) {
            if(t==1)
			cutOpt.add(TC[t][k] >= S[t][0][k]*Y[t][0][k] + tpcx[t][k] * x[t][k]) ;
		  }
        }




		/*****************************************************************************************/
		///           CONTA O TEMPO DE CORTE DO PRIMEIRO PADRAO DE CORTE DO PERIODO  T!=1     /////
		/*****************************************************************************************/

/*		for (int t = 1; t <= periodos; t++) {
			for (int k = 1; k<=numMaxPadroes; k++) {

            IloExpr somatot(env);
             for (int j = 1; j<=numMaxPadroes; j++) {
                somatot += S[t][j][k]*Y[t][j][k];
             }

            cutOpt.add(TC[t][k] >=  tpcx[t][k] * x[t][k] + somatot);
			}
		}


		/*****************************************************************************************/
		///   CONTA O TEMPO DE CORTE " A PARTIR " DO PRIMEIRO PADRAO DE CORTE DO PERIODO      /////
       /*****************************************************************************************/

/*		for (int t = 1; t <= periodos; t++) {
			for (int j = 1; j <=numMaxPadroes; j++) {
				for (int k = 1; k <=numMaxPadroes; k++) {
					if(k!=j)
						cutOpt.add(TC[t][k] >= TC[t][j] + S[t][j][k] + tpcx[t][k]* x[t][k]  - M*( 1 - Y[t][j][k]));
				}
			}
		}




		for (int t = 1; t <= periodos; t++) {
			for (int j = 1; j <=numMaxPadroes; j++) {
				cutOpt.add(TP[t] >= TC[t][j]);
			}
		}




	 	for (int t = 1; t <= periodos; t++) {
			IloExpr soma1(env);
			IloExpr soma2(env);

        for (int k = 1; k <=numMaxPadroes; k++) {
            soma1+= tpcx[t][k]* x[t][k];


            if(t==1)
            soma2 += S[t][0][k] * Y[t][0][k];


           for (int j = 1; j <=numMaxPadroes; j++) {
             if(k!=j)
                soma2 +=  S[t][j][k] * Y[t][j][k];
				}
			}
			cutOpt.add(TP[t] == soma1 + soma2);
		}




		/**** Restrição para garantir que o Tempo de corte só é contabilizado para padrões cortados */
		///Nova
	/*	for (int t = 1; t <= periodos; t++) {
			for (int j = 1; j <=numMaxPadroes; j++) {
				cutOpt.add(TC[t][j] <= M * v[t][j]);
			}
		}



		for (int t = 1; t <= periodos; t++) {
			cutOpt.add(TC[t][0] == 0);
		}






/**************************************************************************************************/
///     ADICIONA COLUNAS INICIAIS : COLUNAS HOMOGENEAS                                           ///
/**************************************************************************************************/


   for (IloInt t = 1; t <= periodos ; t++){
        for (IloInt i= 1; i<=numItens; i++){

          IloNumColumn col = RollsUsed(1);

       for (IloInt j= 1; j<=numItens; j++){
                col += Fill[t][i](PadraoCorte1[j][i]);
                col += FillSetup[t][i](1);
       x[t].add(IloNumVar(col, 0, IloInfinity));


        IloNumColumn col2 = RollsUsed(0);
                col2 += FillSetup[t][i](-M);
          v[t].add(IloNumVar(col2, 0, 1));

            }
       }

   }



        IloCplex cutSolver(cutOpt);

        cutSolver.setParam( IloCplex::Param::MIP::Tolerances::MIPGap, 10e-12);
        cutSolver.setOut(env.getNullStream());
        cutSolver.setWarning(env.getNullStream());



        IloModel patGen (env);
   ///----------------------------------------///
   ///        DADOS PARA A MOCHILA            ///
   ///----------------------------------------///

        IloArray<IloNumVarArray> Use(env, periodos+1);
        for(IloInt t=1; t<=periodos; t++)
            Use[t] = IloNumVarArray(env, numItens+1, 0, IloInfinity, ILOINT);


        IloArray<IloNumArray> newPatt(env, periodos+1);
        for(IloInt t=1; t<=periodos; t++)
            newPatt[t] = IloNumArray(env, numItens+1);

        IloArray<IloNumArray> priceCortes(env, periodos+1);
        for(IloInt t=1; t<=periodos; t++)
            priceCortes[t] = IloNumArray(env, numItens+1);




        IloObjective ReducedCost = IloAdd(patGen, IloMinimize(env));


        for(IloInt t=1; t<=periodos; t++)
         patGen.add(IloScalProd(size, Use[t]) <= CompObjeto);

        IloCplex patSolver(patGen);

        patSolver.setParam( IloCplex::Param::MIP::Tolerances::MIPGap, 10e-12);
        patSolver.setOut(env.getNullStream());
        patSolver.setWarning(env.getNullStream());


        /**************************************************************/
        /********        PROCESSO DE GERAÇÃO DE COLUNAS      **********/
        /**************************************************************/

/*
for(;;){


            cutSolver.solve();

            report1 (cutSolver, periodos, x,   Fill);


            IloNumArray novaColuna(env,periodos+1, 0, IloInfinity, ILOINT);
            IloInt NumNovaCol = 0; // Contador para identificar quantas novas colunas serão adicionadas

      //  for(IloInt t=1; t<=periodos; t++){
      //      novaColuna[t] = 0; // Flag para identificar se haverá uma coluna para o período t ou não
       // }



           for(IloInt t=1; t<=periodos; t++){
                for(IloInt i=1; i<=numItens ; i++)
                    priceCortes[t][i] = - cutSolver.getDual(Fill[t][i]);
            }

            for(IloInt t=1; t<=periodos; t++){
                IloExpr objExpr = ReducedCost.getExpr();///
                objExpr.clear();///

               objExpr += 1 ;

                for(IloInt i=1 ; i<=numItens; i++)
                    objExpr += Use[t][i]*priceCortes[t][i];



                ReducedCost.setExpr( objExpr);//
                objExpr.end();

                patSolver.solve();

                report2 (patSolver,periodos, Use, ReducedCost);
            }

    //cout << " Nova Coluna = " << novaColuna[t] << endl;


/*
              if (patSolver.getValue(ReducedCost) > - RC_EPS){

                        cout <<  " ---- Nesse periodo nao existe padrao que melhore! ---- " << endl;
                    	cout << " PARE! Custo Reduzido = " << patSolver.getValue(ReducedCost) << endl;
                        cout << " Nova Coluna = " << novaColuna[t] << endl;
                        cout << " "  << endl;

                    continue;  //Elisama Comentario: esse parte eh necessaria, pois coso contrario se add uma coluna que ja foi add
                }
                else {
                novaColuna[t] = novaColuna[t] + 1 ; // Contabiliza como 1 caso uma coluna para t seja zerada
               // NumNovaCol = NumNovaCol + 1; // incrementa o número de colunas que serão adicionadas
                }


             //   patSolver.getValues(newPatt[t], Use[t]);
               cout  << " padrão de corte = " << newPatt[t] << " -- j -- "<< x[t].getSize() + 1 << endl;
               // cout << " Número de Novas Colunas :: " << NumNovaCol << endl;
*/
      //    }

 //  if(NumNovaCol == 0){  // Para a execução apenas se o número de novas colunas for zero
 //     break;
  // }


/*
    for(IloInt t=1; t<=periodos; t++){


//       if(novaColuna[t] == 1){  //Aqui: só adiciona uma coluna para o período t se o custo reduzido foi negativo, ou seja, chegou no else da linha 521


         IloNumColumn col = RollsUsed(1);


            for(IloInt i=1; i<=numItens; i++){
                col += Fill[t][i](newPatt[t][i]);
            }

                col += FillSetup[t][x[t].getSize()](1);

        x[t].add(IloNumVar(col, 0, IloInfinity));



         IloNumColumn col2 = RollsUsed(0);
                col2 += FillSetup[t][x[t].getSize()-1](-M);
         v[t].add(IloNumVar(col2, 0, 1));


       //  } // if
      }  // foR
 */
 // } // {;;}


/*
       for(IloInt t=1; t<=periodos; t++){
                   for(IloInt j=1; j<= x[t].getSize(); j++){
                    cutOpt.add(IloConversion(env, x[t][j], ILOINT));
                    cutOpt.add(IloConversion(env, v[t][j], ILOBOOL));
                   }
       }
*/

        cutSolver.setParam(IloCplex::Threads,1);
        cutSolver.setParam(IloCplex::Param::TimeLimit,1200);


        cutSolver.solve();
        cutSolver.exportModel("Model-ModelColunas2.lp");

        cutSolver.setOut(env.getNullStream());
        cutSolver.setWarning(env.getNullStream());
/*
        IloExpr NumPadroes(env);
        IloExpr setups(env);
        IloExpr objObjetosCortados(env);
        IloExpr objFOTRgc(env);
        IloExpr objFOARgc(env);
        IloExpr objFOgc(env);
        IloExpr gap_CPLEX(env);



        for(IloInt t=1; t<= periodos; t++){
            NumPadroes += x[t].getSize() - numItens;

            for(IloInt j=1; j<= x[t].getSize(); j++){
                objObjetosCortados += cutSolver.getValue(x[t][j]);
                setups += cutSolver.getValue(v[t][j]);
            }

            objFOTRgc += cutSolver.getValue(T[t]);
            objFOARgc += cutSolver.getValue(E[t]);
//            objFOgc += cutSolver.getValue(TempoPreparo[t]);
            // gap_CPLEX += cutSolver.
        }

        IloNumExprArg MediaNumPadroes = NumPadroes / periodos;
*/
        IloNum OBJFOGC = cutSolver.getObjValue();



                  crono.stop();
                  double time_crono = crono.getTime() ;

      cout << " função " << OBJFOGC << endl;
      cout << " time "  << time_crono << endl;


      report3 (cutSolver,
		 periodos,
		numItens,
         S,
		 T,
		E,
		TP,
		 TI,
        EstoqueItens,
		 x,
		v,
		 Y);


}



    catch (IloException& ex)
    {
        cerr << "Erro Cplex: " << ex << endl;
    }

    catch (...)
    {
        cerr << "Erro Cpp:" << endl;
    }


    env.end();
}



int main(){

    char instancia[100];

    for(int i=1; i<=1; i++) {

       sprintf(instancia, "/home/prof/Documentos/Elisama-2024/heiristica-ffd-Melega/Testinho.txt ", i );


        principal(instancia);
    }

    return 0;
}


static void report1 (IloCplex& cutSolver,
                     IloNum& periodos,
                     IloArray<IloNumVarArray> Cut,
                     IloArray<IloRangeArray> Fill)
{

}


static void report2 (IloAlgorithm& patSolver,
                     IloNum& periodos,
                     IloArray<IloNumVarArray> Use,
                     IloObjective obj)
{

for(IloInt t=1; t<=periodos; t++){
 cout << endl;
   cout << "Reduced cost is " << patSolver.getValue(obj) << endl;
   cout << endl;
   if (patSolver.getValue(obj) <= -RC_EPS) {
      for (IloInt i = 0; i < Use[t].getSize(); i++)  {
         cout << "  Use" << i << " = " << patSolver.getValue(Use[t][i]) << endl;
      }
      cout << endl;
   }
 }
}


/** apresenta o resultado da funcao objetivo  **/
static void report3 (IloCplex& cutSolver,
		IloNum& periodos,
		IloNum& numItens,
        IloArray<IloNumMatrix2> S,
		IloNumVarArray T,
		IloNumVarArray E,
		IloNumVarArray TP,
		IloNumVarArray TI,
        IloArray<IloNumVarArray>EstoqueItens,
		IloArray<IloNumVarArray> x,
		IloArray<IloNumVarArray> v,
		IloArray<NumVarMatrix2> Y){
}



