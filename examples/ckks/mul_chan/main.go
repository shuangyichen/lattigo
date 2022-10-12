package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func ReadCsv(filename string, num int) (x [][]complex128, y []int) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	x = make([][]complex128, num)
	y = make([]int, num)
	fmt.Println("Sample number", num)

	for i := 0; i < num; i++ {
		x[i] = make([]complex128, 784)
		for j := 0; j < 785; j++ {
			if j == 0 {
				value, err := strconv.Atoi(Rec[i][j])
				if err != nil {
					fmt.Println(err)
				}
				y[i] = value
			} else {
				value, err := strconv.Atoi(Rec[i][j])
				if err != nil {
					fmt.Println(err)
				}
				x[i][j-1] = complex(float64(value)/float64(255), 0)

			}
		}
	}
	return x, y
}
func Readcsv_filter(filename string, f_h int, f_w int, c_in int, c_out int) (filters [][][]complex128) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	filters = make([][][]complex128, c_out)
	for i := 0; i < c_out; i++ {
		filters[i] = make([][]complex128, c_in)
		for j := 0; j < c_in; j++ {
			filters[i][j] = make([]complex128, f_h*f_w)
		}
	}
	for i := 0; i < len(Rec); i++ {
		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}

		x := i / (c_in * c_out) //filter_idx
		tmp := i % (c_in * c_out)
		y := tmp / c_out
		z := tmp % c_out
		// z := i / (c_in * f_h * f_w)
		// y := (i - z*(c_in*f_h*f_w)) / (f_h * f_w)
		// x := i % (f_h * f_w)

		filters[z][y][x] = complex(value, 0)
	}
	return filters
}

func Readcsv_bias(filename string, c_in int, c_out int) (bias []float64) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	bias = make([]float64, c_out)

	for i := 0; i < len(Rec); i++ {
		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}

		bias[i] = value //complex(value, 0)
	}
	return bias
}

func Readcsv_fc_bias(filename string, output int) (bias []float64) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	bias = make([]float64, output)
	for i := 0; i < len(Rec); i++ {
		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}

		bias[i] = value //complex(value, 0)
	}
	return bias
}

func Readcsv_fc(filename string, output int, features_num int) (fc [][]complex128) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	fc = make([][]complex128, output)

	for i := 0; i < output; i++ {
		fc[i] = make([]complex128, features_num)
	}

	for i := 0; i < len(Rec); i++ {
		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}
		j := i / output
		m := i % output

		fc[m][j] = complex(value, 0)
		// fmt.Println(fc[m][j])
	}
	// fmt.Println(fc[:][0])
	return fc
}
func encoding_fc(output int, channel_in int, H int, W int, size int, fea_num int, fc_value [][]complex128, encoder ckks.Encoder, params ckks.Parameters) (fc_pt [][]*ckks.Plaintext) {
	fc_mat := make([][][]complex128, output)
	for i := 0; i < output; i++ {
		fc_mat[i] = make([][]complex128, channel_in)
		for j := 0; j < channel_in; j++ {
			fc_mat[i][j] = make([]complex128, size)
		}
	}
	// fmt.Println(channel_in)
	// for i := 0; i < output; i++ {
	// 	for j := 0; j < channel_in; j++ {
	// 		for k := 0; k < W; k++ {
	// 			for m := 0; m < H; m++ {
	// 				// fmt.Println(j)
	// 				fc_mat[i][j][k*36*5+m*36] = fc_value[i][j*W+k*H+m]
	// 			}
	// 		}
	// 	}
	// }
	// for i := 0; i < output; i++ {
	// 	for j := 0; j < fea_num; j++ {
	// 		k := j / 25
	// 		l := j - 25*k
	// 		m := l / 5
	// 		n := l % 5
	// 		fc_mat[i][k][(m*13+n*2)*18] = fc_value[i][j]

	// 	}
	// }
	for k := 0; k < output; k++ { //10
		for q := 0; q < channel_in; q++ { //4
			for i := 0; i < 20; i++ {
				if i%4 == 0 {
					for j := 0; j < 5; j++ {
						// fmt.Println()
						fc_mat[k][q][i*13*18+j*36] = fc_value[k][4*((i/4)*5+j)+q]
						// fmt.Println(fc_mat[k][q][i*13*18+j*36])
					}
				}
			}
		}
	}
	fc_pt = make([][]*ckks.Plaintext, output)
	for m := 0; m < output; m++ {
		fc_pt[m] = make([]*ckks.Plaintext, channel_in)
		for j := 0; j < channel_in; j++ {
			fc_pt[m][j] = encoder.EncodeNew(fc_mat[m][j], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}

	return fc_pt
}

func encoding_filter_bias(c_in int, c_out int, filters int, filter_size int, H_out int, W_out int, H_in int, W_in int, total_size int, filter_value [][][]complex128, encoder ckks.Encoder, params ckks.Parameters) (pt_w [][]*ckks.Plaintext) {
	enc_filter := make([][][]complex128, c_out)
	// enc_bias := make([][]complex128, c_out)
	if filters == 9 {
		for m := 0; m < c_out; m++ {
			enc_filter[m] = make([][]complex128, c_in)
			// enc_bias[m] = make([][]complex128, c_in)
			for l := 0; l < c_in; l++ {
				enc_filter[m][l] = make([]complex128, total_size)
				// enc_bias[m][l] = make([]complex128, H_out*W_out*filters)
				for i := 0; i < H_out; i++ {
					for j := 0; j < W_out; j++ {
						for k := 0; k < filter_size; k++ {
							for n := 0; n < filter_size; n++ {
								enc_filter[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = filter_value[m][l][k*3+n]
								// if k == 0 && n == 0 {
								// 	enc_bias[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = bias_value[m][l]
								// } else {
								// 	enc_bias[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
								// }
							}
						}
					}
				}
			}
		}
	} else if filters == 18 {
		// fmt.Println(H_out)
		// fmt.Println(W_out)
		// fmt.Println(filters)
		for m := 0; m < c_out; m++ {
			enc_filter[m] = make([][]complex128, c_in)
			// enc_bias[m] = make([][]complex128, c_in)
			for l := 0; l < c_in; l++ {
				enc_filter[m][l] = make([]complex128, total_size)
				// enc_bias[m][l] = make([]complex128, H_out*W_out*filters)
				for i := 0; i < H_out; i++ {
					for j := 0; j < W_out; j++ {
						for k := 0; k < filter_size; k++ {
							for n := 0; n < filter_size; n++ {
								// if i%2 == 0 {
								// fmt.Println(filter_value[m][l][k*3+n])
								enc_filter[m][l][(H_in*(i)*2+j)*filters+k*filter_size+n] = filter_value[m][l][k*3+n]
								// }

							}
						}
					}
				}
			}
		}
		// fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f  %6.10f %6.10f  %6.10f  ...\n", enc_filter[0][0][0], enc_filter[0][0][1], enc_filter[0][0][2], enc_filter[0][0][8], enc_filter[0][0][9], enc_filter[0][0][10], enc_filter[0][0][18])

	}

	pt_w = make([][]*ckks.Plaintext, c_out)
	// pt_b = make([][]*ckks.Plaintext, c_out)
	for m := 0; m < c_out; m++ {
		// fmt.Println(m)
		pt_w[m] = make([]*ckks.Plaintext, c_in)
		// pt_b[m] = make([]*ckks.Plaintext, c_in)
		for l := 0; l < c_in; l++ {
			// fmt.Println(l)
			pt_w[m][l] = encoder.EncodeNew(enc_filter[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			// pt_b[m][l] = encoder.EncodeNew(enc_bias[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}
	return pt_w
}

func main() {

	flag.Parse()

	var err error

	// var btp *bootstrapping.Bootstrapper
	var kgen rlwe.KeyGenerator
	var encoder ckks.Encoder
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var encryptor ckks.Encryptor
	var decryptor ckks.Decryptor
	// var plaintext *ckks.Plaintext

	// r, _ := gonpy.NewFileReader("00.npy")
	// fmt.Println(r.Shape)
	// data, _ := r.GetFloat64()

	// for i := 0; i < r.Shape[0]; i++ {
	// 	fmt.Println(data)
	// }

	// Bootstrapping parameters
	// Two sets of four parameters each, DefaultParametersSparse and DefaultParametersDense,
	// (each index 0 to 3) ensuring 128 bit of security are available in
	// github.com/tuneinsight/lattigo/v3/ckks/bootstrapping/default_params.go
	//
	// LogSlots is hardcoded to 15 in the parameters, but can be changed from 1 to 15.
	// When changing LogSlots make sure that the number of levels allocated to CtS and StC is
	// smaller or equal to LogSlots.

	F_H := 28
	F_W := 28
	party1_cols := 1
	//vector<int> party1_index(F_H_out);
	//vector<int> party2_index(F_H*F_W);

	filter_size := 3
	stride := 1
	F_H_out := (F_H - filter_size + 1) / stride
	F_W_out := (F_W - filter_size + 1) / stride
	F_H_out_1 := 11
	F_W_out_1 := 11
	filters := filter_size * filter_size
	paramSet := bootstrapping.DefaultParametersSparse[3] // bootstrapping.DefaultParametersDense[0]
	ckksParams := paramSet.SchemeParams

	if *flagShort {
		ckksParams.LogN = 13
		ckksParams.LogSlots = 12
	}

	btpParams := paramSet.BootstrappingParams

	params, err := ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, H(%d; %d), logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), params.HammingWeight(), btpParams.EphemeralSecretWeight, params.LogQP(), params.QCount(), math.Log2(params.DefaultScale()), params.Sigma())
	prediction_correct := 0
	rots := []int{1, 2, 3, 4, 5, 6, 7, 8, F_H_out * F_W_out * filters, 9, 234, 243, 17, 34, 465, 482, 499, 930, 947, 964, 198, 216, 36, 234, 252, 18, 468, 486, 936}
	// rots_
	// rots := make([]int, filters+1)
	// for i := 0; i < filters; i++ {
	// 	rots[i] = i
	// }
	// rots[filters] = F_H_out * F_W_out * filters

	kgen = ckks.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPair()
	rlk := kgen.GenRelinearizationKey(sk, 2)
	rotk := kgen.GenRotationKeysForRotations(append(rots, 0-filters), true, sk)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotk})

	encoder = ckks.NewEncoder(params)
	decryptor = ckks.NewDecryptor(params, sk)
	encryptor = ckks.NewEncryptor(params, pk)

	c_in_0 := 1
	c_out_0 := 4
	c_in_1 := 4
	c_out_1 := 4
	features_num := 100
	// c_in_1 := 32
	// c_out_1 := 64
	filters0_value := Readcsv_filter("00.csv", filter_size, filter_size, c_in_0, c_out_0)
	// fmt.Println(filters0_value[0][0])
	// fmt.Println(filters0_value[1][0])
	// fmt.Println(filters0_value[2][0])
	// fmt.Println(filters0_value[3][0])

	bias0_value := Readcsv_bias("01.csv", c_in_0, c_out_0)
	filters1_value := Readcsv_filter("30.csv", filter_size, filter_size, c_in_1, c_out_1)
	bias1_value := Readcsv_bias("31.csv", c_in_1, c_out_1)
	fc_value := Readcsv_fc("70.csv", 10, features_num)

	fc_bias := Readcsv_fc_bias("71.csv", 10)
	pt_w0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, F_H_out+2, F_W_out+2, F_H_out*F_W_out*filters, filters0_value, encoder, params)

	// fmt.Printf("Encoding finished")
	filters_1 := 18
	pt_w1 := encoding_filter_bias(c_in_1, c_out_1, filters_1, filter_size, F_H_out_1, F_W_out_1, F_H_out_1+2, F_W_out_1+2, F_H_out*F_W_out*filters, filters1_value, encoder, params)
	pt_fc := encoding_fc(10, 4, 5, 5, 26*26*9, 100, fc_value, encoder, params)
	fmt.Printf("Encoding finished \n")

	// Scheme context and keys
	//Preparing data
	sample_nums := 1000
	sample_num := 1
	xs, ys := ReadCsv("mnist_test.csv", sample_nums)
	for i := 0; i < 100; i++ {

		// x := xs[1]
		x := make([][]complex128, sample_num)
		x[0] = xs[400+i]
		label := ys[400+i]
		// x[0] = make([]complex128, 784)
		// for i := 0; i < len(x[0]); i++ {
		// 	x[0][i] = complex(float64(i), 0)
		// }
		// fmt.Println(len(x))
		// fmt.Println(len(x[0]))
		// // fmt.Println(len(y))

		// F_H := 28
		// F_W := 28
		// party1_cols := 1
		// //vector<int> party1_index(F_H_out);
		// //vector<int> party2_index(F_H*F_W);

		// filter_size := 3
		// stride := 1
		// F_H_out := (F_H - filter_size + 1) / stride
		// F_W_out := (F_W - filter_size + 1) / stride
		// F_H_out_1 := 11
		// F_W_out_1 := 11
		// filters := filter_size * filter_size

		features1 := make([][]complex128, sample_num)
		features2 := make([][]complex128, sample_num)
		for i := 0; i < sample_num; i++ {
			features1[i] = make([]complex128, F_H_out*F_W_out*filters)
			features2[i] = make([]complex128, F_H_out*F_W_out*filters)
		}

		//Setting filter 0
		// c_in_0 := 1
		// c_out_0 := 4
		// c_in_1 := 4
		// c_out_1 := 4
		// features_num := 100
		// // c_in_1 := 32
		// // c_out_1 := 64
		// filters0_value := Readcsv_filter("00.csv", filter_size, filter_size, c_in_0, c_out_0)
		// // fmt.Println(filters0_value[0][0])
		// // fmt.Println(filters0_value[1][0])
		// // fmt.Println(filters0_value[2][0])
		// // fmt.Println(filters0_value[3][0])

		// bias0_value := Readcsv_bias("01.csv", c_in_0, c_out_0)
		// filters1_value := Readcsv_filter("30.csv", filter_size, filter_size, c_in_1, c_out_1)
		// bias1_value := Readcsv_bias("31.csv", c_in_1, c_out_1)
		// fc_value := Readcsv_fc("70.csv", 10, features_num)

		// fc_bias := Readcsv_fc_bias("71.csv", 10)
		// pt_w0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, F_H_out+2, F_W_out+2, F_H_out*F_W_out*filters, filters0_value, encoder, params)

		// // fmt.Printf("Encoding finished")
		// filters_1 := 18
		// pt_w1 := encoding_filter_bias(c_in_1, c_out_1, filters_1, filter_size, F_H_out_1, F_W_out_1, F_H_out_1+2, F_W_out_1+2, F_H_out*F_W_out*filters, filters1_value, encoder, params)
		// pt_fc := encoding_fc(10, 4, 5, 5, 26*26*9, 100, fc_value, encoder, params)
		// fmt.Printf("Encoding finished \n")

		// pt_w0, pt_b0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, filters0_value, bias0_value, encoder, params)
		// pt_w1, pt_b1 := encoding_filter_bias(c_in_1, c_out_1, filters, filter_size, F_H_out_1, F_W_out_1, filters1_value, bias1_value, encoder, params)
		// fmt.Printf("Read filter and bias \n")
		// filters0_value := make([][][]complex128, c_out_0)
		// bias0_value := make([][]complex128, c_out_0)
		// for i := 0; i < c_out_0; i++ {
		// 	filters0_value[i] = make([][]complex128, c_in_0)
		// 	bias0_value[i] = make([]complex128, c_in_0)
		// 	for j := 0; j < c_in_0; j++ {
		// 		filters0_value[i][j] = make([]complex128, filters)
		// 		bias0_value[i][j] = complex(1.0, 0)
		// 		for p := 0; p < filters; p++ {
		// 			filters0_value[i][j][p] = complex(1.0, 0)
		// 		}
		// 	}
		// }

		for i := 0; i < F_H_out; i++ {
			for j := 0; j < F_W_out; j++ {
				for k := 0; k < filter_size; k++ {
					for n := 0; n < filter_size; n++ {
						if j+n < party1_cols {
							features1[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = x[0][(i+k)*F_H+j+n]
							features2[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
						} else {
							features2[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = x[0][(i+k)*F_H+j+n]
							features1[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
						}
					}
				}
			}
		}
		// for i := 0; i < 26*26*9; i++ {
		// 	if real(features2[0][i]) > 0.7 {
		// 		fmt.Println(i)
		// 		fmt.Println(real(features2[0][i]))
		// 	}
		// }
		// for i := 0; i < 28*28; i++ {
		// 	if real(x[0][i]) > 0.7 {
		// 		fmt.Println(i)
		// 		fmt.Println(real(x[0][i]))
		// 	}
		// }
		// v := float64(76) / float64(255)
		// fmt.Println(v)

		//generating filter
		// for i := 0; i < filters; i++ {
		// 	// filters_value[i] = complex(rand.Float64(), 0)
		// 	filters_value[i] = complex(1.0, 0)
		// }

		//encoding filter 0
		// enc_filter0 := make([][][]complex128, c_out_0)
		// enc_bias0 := make([][][]complex128, c_out_0)
		// for m := 0; m < c_out_0; m++ {
		// 	enc_filter0[m] = make([][]complex128, c_in_0)
		// 	enc_bias0[m] = make([][]complex128, c_in_0)
		// 	for l := 0; l < c_in_0; l++ {
		// 		enc_filter0[m][l] = make([]complex128, F_H_out*F_W_out*filters)
		// 		enc_bias0[m][l] = make([]complex128, F_H_out*F_W_out*filters)
		// 		for i := 0; i < F_H_out; i++ {
		// 			for j := 0; j < F_W_out; j++ {
		// 				for k := 0; k < filter_size; k++ {
		// 					for n := 0; n < filter_size; n++ {
		// 						enc_filter0[m][l][(F_W_out*(i)+j)*filters+k*filter_size+n] = filters0_value[m][l][k*3+n]
		// 						if k == 0 && n == 0 {
		// 							enc_filter0[m][l][(F_W_out*(i)+j)*filters+k*filter_size+n] = bias0_value[m][l]
		// 						} else {
		// 							enc_filter0[m][l][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }

		// rots := []int{1, 2, 3, 4, 5, 6, 7, 8, F_H_out * F_W_out * filters, 9, 234, 243, 17, 34, 465, 482, 499, 930, 947, 964, 198, 216, 36, 234, 252, 18, 468, 486, 936}
		// // rots_
		// // rots := make([]int, filters+1)
		// // for i := 0; i < filters; i++ {
		// // 	rots[i] = i
		// // }
		// // rots[filters] = F_H_out * F_W_out * filters

		// kgen = ckks.NewKeyGenerator(params)

		// sk, pk = kgen.GenKeyPair()
		// rlk := kgen.GenRelinearizationKey(sk, 2)
		// rotk := kgen.GenRotationKeysForRotations(append(rots, 0-filters), true, sk)
		// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotk})

		// encoder = ckks.NewEncoder(params)
		// decryptor = ckks.NewDecryptor(params, sk)
		// encryptor = ckks.NewEncryptor(params, pk)

		// pt_w0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, F_H_out+2, F_W_out+2, F_H_out*F_W_out*filters, filters0_value, encoder, params)

		// // fmt.Printf("Encoding finished")
		// filters_1 := 18
		// pt_w1 := encoding_filter_bias(c_in_1, c_out_1, filters_1, filter_size, F_H_out_1, F_W_out_1, F_H_out_1+2, F_W_out_1+2, F_H_out*F_W_out*filters, filters1_value, encoder, params)
		// pt_fc := encoding_fc(10, 4, 5, 5, 26*26*9, 100, fc_value, encoder, params)
		// fmt.Printf("Encoding finished \n")
		//encoding weights

		// pt_w0 := make([][]*ckks.Plaintext, c_out_0)
		// pt_b0 := make([][]*ckks.Plaintext, c_out_0)
		// for m := 0; m < c_out_0; m++ {
		// 	pt_w0[m] = make([]*ckks.Plaintext, c_in_0)
		// 	pt_b0[m] = make([]*ckks.Plaintext, c_in_0)
		// 	for l := 0; l < c_in_0; l++ {
		// 		pt_w0[m][l] = encoder.EncodeNew(enc_filter0[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// 		pt_b0[m][l] = encoder.EncodeNew(enc_filter0[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// 	}
		// }

		// pt_w := encoder.EncodeNew(enc_filter, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

		//encode features
		pt_f1 := make([]*ckks.Plaintext, sample_num)
		pt_f2 := make([]*ckks.Plaintext, sample_num)
		ct_f1 := make([]*ckks.Ciphertext, sample_num)
		ct_f2 := make([]*ckks.Ciphertext, sample_num)
		for i := 0; i < sample_num; i++ {
			pt_f1[i] = encoder.EncodeNew(features1[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			pt_f2[i] = encoder.EncodeNew(features2[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			ct_f1[i] = encryptor.EncryptNew(pt_f1[i])
			ct_f2[i] = encryptor.EncryptNew(pt_f2[i])
		}
		//aggregating clients features

		rots_for_replica_1 := []int{17, 34, 465, 482, 499, 930, 947, 964}
		rots_for_pool_1 := []int{9, 234, 243}
		rots_for_pool_2 := []int{18, 468, 486}
		// rots_for_pool_1 := []int{17, 34, 231, 248, 265, 462, 479, 496}
		ct_result := make([][]*ckks.Ciphertext, sample_num)
		for i := 0; i < sample_num; i++ {
			// printDebug(params, ct_f1[i], decryptor, encoder)
			// printDebug(params, ct_f2[i], decryptor, encoder)
			evaluator.Add(ct_f1[i], ct_f2[i], ct_f1[i])
			// printDebug(params, ct_f1[i], decryptor, encoder)
			t1 := time.Now()
			ct_result[i] = Conv1(params, ct_f1, pt_w0, bias0_value, c_in_0, c_out_0, filters, F_H_out*F_W_out*filters, evaluator, decryptor, encoder, 9, 26*26*9)
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)

			t2 := time.Now()
			fmt.Println("Conv time: ", t2.Sub(t1))
			t3 := time.Now()
			for j := 0; j < c_out_0; j++ {
				ct_result[i][j] = relu(params, ct_result[i][j], evaluator, decryptor, encoder)
			}
			t4 := time.Now()
			fmt.Println("ReLU time: ", t4.Sub(t3))
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)

			// // fmt.Println(ct_result[i].Level())
			// // ct_result[i] = relu(params, ct_result[i], evaluator, decryptor, encoder)
			// t3 := time.Now()
			// fmt.Println("ReLU time: ", t3.Sub(t2))
			// // fmt.Println(ct_result[i].Level())
			t5 := time.Now()
			for j := 0; j < c_out_0; j++ {
				ct_result[i][j] = avg_pooling(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder, rots_for_pool_1, 18, 13*13*18)
			}
			t6 := time.Now()
			fmt.Println("AvgPooling time: ", t6.Sub(t5))
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)

			// t4 := time.Now()
			// // fmt.Println(ct_result[i].Level())
			// fmt.Println("Pooling time: ", t4.Sub(t3))
			// // printDebug(params, ct_result[i], decryptor, encoder)
			t7 := time.Now()
			for j := 0; j < c_out_0; j++ {
				ct_result[i][j] = rot_for_rep_conv(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder, rots_for_replica_1)
			}
			t8 := time.Now()
			fmt.Println("Rot for replica time: ", t8.Sub(t7))
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)
			// t5 := time.Now()
			// fmt.Println("Replication time: ", t5.Sub(t4))
			// ct_result[i].Rescale()
			// printDebug(params, ct_result[i], decryptor, encoder)
			//conv again
			fmt.Println("**************Second layer**********************")
			t9 := time.Now()
			ct_result[i] = Conv2(params, ct_result[i], pt_w1, bias1_value, c_in_1, c_out_1, filters, F_H_out*F_W_out*filters, evaluator, decryptor, encoder, 18, 13*13*18)
			t10 := time.Now()
			fmt.Println("Second conv time: ", t10.Sub(t9))
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)
			t11 := time.Now()
			for j := 0; j < c_out_1; j++ {
				ct_result[i][j] = relu(params, ct_result[i][j], evaluator, decryptor, encoder)
			}
			t12 := time.Now()
			fmt.Println("ReLU time: ", t12.Sub(t11))
			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)
			t13 := time.Now()
			for j := 0; j < c_out_1; j++ {
				ct_result[i][j] = avg_pooling2(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder, rots_for_pool_2, 18, 13*13*18)
			}
			t14 := time.Now()
			fmt.Println("AvgPooling time: ", t14.Sub(t13))
			// printDebug2(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// printDebug(params, ct_result[i][3], decryptor, encoder)
			t15 := time.Now()
			result := FC(params, ct_result[i], 10, pt_fc, fc_bias, 4, evaluator, decryptor, encoder)
			t16 := time.Now()
			fmt.Println("FC time: ", t16.Sub(t15))

			prediction_res := Decrypt(params, result, decryptor, encoder)
			fmt.Println("Prediction result: ", prediction_res)
			fmt.Println("Label: ", label)
			if prediction_res == label {
				prediction_correct += 1
				// fmt.Println("Prediction result: ", prediction_res)
				// fmt.Println("Label: ", label)
				fmt.Println("Correct Prediction:", prediction_correct)
			}
			// printResult(params, result[0], decryptor, encoder)
			// printResult(params, result[1], decryptor, encoder)
			// printResult(params, result[2], decryptor, encoder)
			// printResult(params, result[3], decryptor, encoder)
			// printResult(params, result[4], decryptor, encoder)
			// printResult(params, result[5], decryptor, encoder)
			// printResult(params, result[6], decryptor, encoder)
			// printResult(params, result[7], decryptor, encoder)
			// printResult(params, result[8], decryptor, encoder)
			// printResult(params, result[9], decryptor, encoder)

			// printDebug(params, ct_result[i][0], decryptor, encoder)
			// printDebug(params, ct_result[i][1], decryptor, encoder)
			// printDebug(params, ct_result[i][2], decryptor, encoder)
			// for j := 0; j < c_out_1; j++ {
			// 	ct_result[i][j] = rot_for_rep_conv(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder)
			// }
			// fmt.Println("channel: ", len(ct_result[i]))
			// printDebug(params, ct_result[i][2], decryptor, encoder)

		}
	}

}

func Decrypt(params ckks.Parameters, res []*ckks.Ciphertext, decryptor ckks.Decryptor, encoder ckks.Encoder) (max int) {
	res_pt := make([]float64, 10)
	for i := 0; i < 10; i++ {
		values := encoder.Decode(decryptor.DecryptNew(res[i]), params.LogSlots())
		res_pt[i] = real(values[0])
		fmt.Println(res_pt[i])
	}
	// fmt.Println()
	max = 0
	for i := 1; i < 10; i++ {
		if res_pt[i] > res_pt[max] {
			max = i
		}
	}
	return max
}
func FC(params ckks.Parameters, data []*ckks.Ciphertext, output int, fc_pt [][]*ckks.Plaintext, bias []float64, in_channel int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result_ct []*ckks.Ciphertext) {
	result_ct = make([]*ckks.Ciphertext, output)
	for i := 0; i < output; i++ {
		result_ct[i] = combine_channels_fc(params, data, fc_pt[i], in_channel, evaluator, decryptor, encoder)
		evaluator.AddConst(result_ct[i], bias[i], result_ct[i])
		// evaluator.Rescale(result_ct[i], params.DefaultScale(), result_ct[i])
	}
	return result_ct
}

func combine_channels_fc(params ckks.Parameters, data []*ckks.Ciphertext, fc_value []*ckks.Plaintext, in_channel int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result *ckks.Ciphertext) {
	result = fc(params, data[0], fc_value[0], 36, 24, evaluator, decryptor, encoder)
	for i := 1; i < in_channel; i++ {
		tmp2 := fc(params, data[i], fc_value[i], 36, 24, evaluator, decryptor, encoder)
		evaluator.Add(result, tmp2, result)
	}
	return result
}

func fc(params ckks.Parameters, data *ckks.Ciphertext, fc_pt *ckks.Plaintext, rot_steps int, rot_times int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (mulres *ckks.Ciphertext) {
	mulres = evaluator.MulNew(data, fc_pt)
	evaluator.Rescale(mulres, params.DefaultScale(), mulres)
	tmp := evaluator.RotateNew(mulres, 36)
	evaluator.Add(mulres, tmp, mulres)
	for i := 1; i < 5; i++ {
		evaluator.Rotate(tmp, 36, tmp)
		evaluator.Add(mulres, tmp, mulres)
	}
	tmp1 := evaluator.RotateNew(mulres, 936)
	evaluator.Add(mulres, tmp1, mulres)
	for i := 2; i < 5; i++ {
		evaluator.Rotate(tmp1, 936, tmp1)
		evaluator.Add(mulres, tmp1, mulres)
	}
	// for i := 0; i < 13; i++ {
	// 	for j := 0; j < 5; j++ {
	// 		mask[i*13*step+j*2*step] = complex(1, 0)
	// 	}
	// }
	return mulres
}

func rot_for_rep_conv(params ckks.Parameters, data *ckks.Ciphertext, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder, rots_for_pool []int) (ct *ckks.Ciphertext) {
	// rots_for_pool := []int{17, 34, 231, 248, 265, 462, 479, 496}
	dup := evaluator.RotateNew(data, size)
	evaluator.Add(dup, data, dup)
	for i := 0; i < len(rots_for_pool); i++ {
		tmp := evaluator.RotateNew(data, rots_for_pool[i])
		evaluator.Add(dup, tmp, dup)
	}

	// mask := make([]complex128, size)
	// for i := 0; i < size; i++ {
	// 	if i%18 == 0 {
	// 		mask[i] = complex(1, 0)
	// 	} else {
	// 		mask[i] = complex(0, 0)
	// 	}
	// }
	// pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	// evaluator.Mul(dup, pt_mask, dup)
	// evaluator.MultByConst(dup, 0.25, dup)
	return dup
}
func avg_pooling(params ckks.Parameters, data *ckks.Ciphertext, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder, rots_for_pool []int, step int, size2 int) (ct *ckks.Ciphertext) {
	// rots_for_pool := []int{9, 234, 243} //for first conv
	dup := evaluator.RotateNew(data, size)
	evaluator.Add(dup, data, dup)
	for i := 0; i < len(rots_for_pool); i++ {
		tmp := evaluator.RotateNew(data, rots_for_pool[i])
		evaluator.Add(dup, tmp, dup)
	}
	mask := make([]complex128, size)
	// for i := 0; i < size2; i++ {
	// 	if i%18 == 0 {
	// 		if (i/(13*18))%2 == 0 {
	// 			mask[i] = complex(1, 0)
	// 		} else {
	// 			mask[i] = complex(0, 0)
	// 		}
	// 	} else {
	// 		mask[i] = complex(0, 0)
	// 	}
	// }
	for i := 0; i < 26; i++ {
		if i%2 == 0 {
			for j := 0; j < 13; j++ {
				mask[i*13*18+j*18] = complex(1, 0)
			}
		}
	}

	pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	evaluator.Mul(dup, pt_mask, dup)
	evaluator.Rescale(dup, params.DefaultScale(), dup)
	evaluator.MultByConst(dup, 0.25, dup)
	evaluator.Rescale(dup, params.DefaultScale(), dup)
	return dup
}

func avg_pooling2(params ckks.Parameters, data *ckks.Ciphertext, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder, rots_for_pool []int, step int, size2 int) (ct *ckks.Ciphertext) {
	// rots_for_pool := []int{9, 234, 243} //for first conv
	dup := evaluator.RotateNew(data, size)
	evaluator.Add(dup, data, dup)
	for i := 0; i < len(rots_for_pool); i++ {
		tmp := evaluator.RotateNew(data, rots_for_pool[i])
		evaluator.Add(dup, tmp, dup)
	}
	mask := make([]complex128, size)

	for i := 0; i < 20; i++ {
		if i%4 == 0 {
			for j := 0; j < 5; j++ {
				mask[i*26*9+j*36] = 1
			}
		}
	}
	// for i := 0; i < 13; i++ {
	// 	for j := 0; j < 5; j++ {
	// 		mask[i*13*step+j*2*step] = complex(1, 0)
	// 	}
	// }
	pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	evaluator.Mul(dup, pt_mask, dup)
	evaluator.Rescale(dup, params.DefaultScale(), dup)
	evaluator.MultByConst(dup, 0.25, dup)
	evaluator.Rescale(dup, params.DefaultScale(), dup)
	return dup
}

func Max(a int, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

func tree_cipher(params ckks.Parameters, data *ckks.Ciphertext, coeffs []float64, degree int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (ct *ckks.Ciphertext) {
	powers := make([]*ckks.Ciphertext, degree)
	// rs := make([]*ckks.Ciphertext, degree+1)
	// fmt.Println(degree)
	powers[0] = data.CopyNew()

	powers[1] = evaluator.MulNew(powers[0], powers[0])
	evaluator.Relinearize(powers[1], powers[1])
	evaluator.Rescale(powers[1], params.DefaultScale(), powers[1])
	// printDebug(params, powers[1], decryptor, encoder)

	for i := 1; i <= degree; i++ {
		// fmt.Println(i)
		evaluator.MultByConst(powers[i-1], coeffs[i], powers[i-1])
		// rs[i] = evaluator.MultByConstNew(powers[i], 1)
		// evaluator.Rescale(powers[i], params.DefaultScale(), powers[i])
	}
	// printDebug(params, powers[1], decryptor, encoder)
	evaluator.Rescale(powers[0], params.DefaultScale(), powers[0])
	evaluator.Rescale(powers[1], params.DefaultScale(), powers[1])
	// evaluator.DropLevel(powers[0], 6)
	// fmt.Println("params")
	// printDebug(params, powers[0], decryptor, encoder)
	// printDebug(params, powers[1], decryptor, encoder)
	ct = evaluator.AddNew(powers[0], powers[1])
	// printDebug(params, ct, decryptor, encoder)
	// for i := 1; i < degree; i++ {
	// 	evaluator.Add(powers[0], powers[i], powers[1])
	// }
	// printDebug(params, ct, decryptor, encoder)
	evaluator.AddConst(ct, coeffs[0], ct)
	// printDebug(params, ct, decryptor, encoder)
	return ct
}
func relu(params ckks.Parameters, data *ckks.Ciphertext, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result *ckks.Ciphertext) {
	degree := 2
	coeffs := make([]float64, degree+1)
	coeffs[0] = 0.00001
	coeffs[1] = 0.00001
	coeffs[2] = 1.0
	result = tree_cipher(params, data, coeffs, degree, evaluator, decryptor, encoder)
	return result
}
func Conv2(params ckks.Parameters, data []*ckks.Ciphertext, weights [][]*ckks.Plaintext, bias []float64, channel_in int, channel_out int, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder, step int, size2 int) (result_ct []*ckks.Ciphertext) {

	result_ct = make([]*ckks.Ciphertext, channel_out)
	for i := 0; i < channel_out; i++ {
		// fmt.Println(i)
		result_ct[i] = combine_channels(params, data, weights[i], channel_in, filters, size, evaluator, decryptor, encoder)
		evaluator.AddConst(result_ct[i], bias[i], result_ct[i])
		mask := make([]complex128, size)
		// for i := 0; i < size2; i++ {
		// 	if i%step == 0 {
		// 		m := (i / step) % 13
		// 		if m == 11 || m == 12 || m == 10 {
		// 			mask[i] = complex(0, 0)
		// 		} else {
		// 			mask[i] = complex(1, 0)
		// 		}
		// 	} else {
		// 		mask[i] = complex(0, 0)
		// 	}
		// }
		for i := 0; i < 21; i++ {
			if i%2 == 0 {
				for j := 0; j < 11; j++ {
					mask[i*26*9+j*18] = 1
				}
			}
		}
		pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		evaluator.Mul(result_ct[i], pt_mask, result_ct[i])
		evaluator.Rescale(result_ct[i], params.DefaultScale(), result_ct[i])
		// evaluator.AddConst(result_ct[i], bias[i], result_ct[i])
	}

	// result_ct[0] = combine_channels(params, data, weights[0], bias[0], channel_in, filters, size, evaluator, decryptor, encoder)
	return result_ct
}

func Conv1(params ckks.Parameters, data []*ckks.Ciphertext, weights [][]*ckks.Plaintext, bias []float64, channel_in int, channel_out int, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder, step int, size2 int) (result_ct []*ckks.Ciphertext) {

	result_ct = make([]*ckks.Ciphertext, channel_out)
	for i := 0; i < channel_out; i++ {
		// fmt.Println(i)
		result_ct[i] = combine_channels(params, data, weights[i], channel_in, filters, size, evaluator, decryptor, encoder)
		evaluator.AddConst(result_ct[i], bias[i], result_ct[i])
		mask := make([]complex128, size)
		for i := 0; i < size2; i++ {
			if i%step == 0 {
				mask[i] = complex(1, 0)
			} else {
				mask[i] = complex(0, 0)
			}
		}
		pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		evaluator.Mul(result_ct[i], pt_mask, result_ct[i])
		evaluator.Rescale(result_ct[i], params.DefaultScale(), result_ct[i])
		// evaluator.AddConst(result_ct[i], bias[i], result_ct[i])
	}

	// result_ct[0] = combine_channels(params, data, weights[0], bias[0], channel_in, filters, size, evaluator, decryptor, encoder)
	return result_ct
}

func combine_channels(params ckks.Parameters, data []*ckks.Ciphertext, weights []*ckks.Plaintext, channel_in int, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result *ckks.Ciphertext) {
	// printDebug(params, data[0], decryptor, encoder)
	// printDebug(params, data[1], decryptor, encoder)
	// printDebug(params, data[2], decryptor, encoder)
	result = conv(params, data[0], weights[0], filters, size, evaluator, decryptor, encoder)
	// printDebug(params, result, decryptor, encoder)
	for j := 1; j < channel_in; j++ {
		tmp := conv(params, data[j], weights[j], filters, size, evaluator, decryptor, encoder)
		// fmt.Println("tmp")
		// fmt.Println(j)
		// printDebug(params, tmp, decryptor, encoder)
		evaluator.Add(result, tmp, result)
		// printDebug(params, result, decryptor, encoder)
	}
	// printDebug(params, result, decryptor, encoder)
	return result
}
func conv(params ckks.Parameters, data *ckks.Ciphertext, weights *ckks.Plaintext, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (dup *ckks.Ciphertext) {
	// printDebug(params, data, decryptor, encoder)
	mulres := evaluator.MulNew(data, weights)
	// printDebug(params, data, decryptor, encoder)
	evaluator.Rescale(mulres, params.DefaultScale(), mulres)
	// printDebug(params, mulres, decryptor, encoder)
	// evaluator.Relinearize(data, data)
	dup = evaluator.RotateNew(mulres, size)
	// values := encoder.Decode(decryptor.DecryptNew(dup), params.LogSlots())
	evaluator.Add(dup, mulres, dup)
	for i := 1; i < filters; i++ {
		tmp := evaluator.RotateNew(mulres, i)
		evaluator.Add(dup, tmp, dup)
	}
	// printDebug(params, dup, decryptor, encoder)
	// mask := make([]complex128, size)
	// for i := 0; i < size; i++ {
	// 	if i%9 == 0 {
	// 		mask[i] = complex(1, 0)
	// 	} else {
	// 		mask[i] = complex(0, 0)
	// 	}
	// }
	// pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	// evaluator.Mul(dup, pt_mask, dup)
	// printDebug(params, dup, decryptor, encoder)
	// evaluator.Relinearize(dup, dup)
	// evaluator.Rescale(dup, params.DefaultScale(), dup)
	// evaluator.Add(dup, bias, dup)
	// printDebug(params, dup, decryptor, encoder)

	// fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", values[0], values[1], values[2], values[3])
	return dup
}

func printDebug(params ckks.Parameters, ciphertext *ckks.Ciphertext, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3], valuesTest[4], valuesTest[5], valuesTest[6], valuesTest[7], valuesTest[8], valuesTest[9], valuesTest[10])
	// fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	// precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	// fmt.Println(precStats.String())
	// fmt.Println(valuesTest[0])
	// fmt.Println(valuesTest[1])
	// fmt.Println(valuesTest[2])
	// fmt.Println(valuesTest[3])
	// fmt.Println(valuesTest[8])
	// fmt.Println(valuesTest[9])
	// fmt.Println(valuesTest[10])
	// fmt.Println(valuesTest[11])
	fmt.Println(valuesTest[0])
	fmt.Println(valuesTest[18])
	fmt.Println(valuesTest[2*13*18])
	fmt.Println(valuesTest[2*13*18+18])
	// fmt.Println(valuesTest[2*26*9])
	// fmt.Println(valuesTest[4*26*9])
	// // fmt.Println(valuesTest[6*26*9])
	// fmt.Println(valuesTest[16*26*9])
	// fmt.Println(valuesTest[20*26*9])
	// fmt.Println(valuesTest[22*26*9])
	// fmt.Println(valuesTest[24*26*9])
	// for i := 0; i < 5000; i++ {
	// 	fmt.Println(valuesTest[i])
	// }

	return
}
func printDebug2(params ckks.Parameters, ciphertext *ckks.Ciphertext, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3], valuesTest[4], valuesTest[5], valuesTest[6], valuesTest[7], valuesTest[8], valuesTest[9], valuesTest[10])
	// fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	// precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	// fmt.Println(precStats.String())
	// fmt.Println(valuesTest[0])
	// fmt.Println(valuesTest[1])
	// fmt.Println(valuesTest[2])
	// fmt.Println(valuesTest[3])
	// fmt.Println(valuesTest[8])
	// fmt.Println(valuesTest[9])
	// fmt.Println(valuesTest[10])
	// fmt.Println(valuesTest[11])
	fmt.Println(valuesTest[0])
	fmt.Println(valuesTest[36])
	fmt.Println(valuesTest[4*13*18])
	fmt.Println(valuesTest[4*13*18+36])
	// fmt.Println(valuesTest[2*26*9])
	// fmt.Println(valuesTest[4*26*9])
	// // fmt.Println(valuesTest[6*26*9])
	// fmt.Println(valuesTest[16*26*9])
	// fmt.Println(valuesTest[20*26*9])
	// fmt.Println(valuesTest[22*26*9])
	// fmt.Println(valuesTest[24*26*9])
	// for i := 0; i < 5000; i++ {
	// 	fmt.Println(valuesTest[i])
	// }

	return
}
func printResult(params ckks.Parameters, ciphertext *ckks.Ciphertext, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale))
	// fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f  %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3], valuesTest[4], valuesTest[5], valuesTest[6], valuesTest[7], valuesTest[8], valuesTest[9], valuesTest[10])
	// fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	// precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	// fmt.Println(precStats.String())
	fmt.Println(valuesTest[0])
	// fmt.Println(valuesTest[18])
	// fmt.Println(valuesTest[36])
	// fmt.Println(valuesTest[198])

	return
}
