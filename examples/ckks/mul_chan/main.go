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

func ReadCsv(filename string, num int) (x [][]complex128, y []complex128) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	x = make([][]complex128, num)
	y = make([]complex128, num)
	fmt.Println("Sample number", num)

	for i := 0; i < num; i++ {
		x[i] = make([]complex128, 784)
		for j := 0; j < 785; j++ {
			if j == 0 {
				value, err := strconv.Atoi(Rec[i][j])
				if err != nil {
					fmt.Println(err)
				}
				y[i] = complex(float64(value), 0)
			} else {
				value, err := strconv.Atoi(Rec[i][j])
				if err != nil {
					fmt.Println(err)
				}
				x[i][0] = complex(float64(value/255), 0)

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
		z := i / (c_in * f_h * f_w)
		y := (i - z*(c_in*f_h*f_w)) / (f_h * f_w)
		x := i % (f_h * f_w)

		filters[z][y][x] = complex(value, 0)
	}
	return filters
}

func Readcsv_bias(filename string, c_in int, c_out int) (bias [][]complex128) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	bias = make([][]complex128, c_out)
	for i := 0; i < c_out; i++ {
		bias[i] = make([]complex128, c_in)
	}
	for i := 0; i < len(Rec); i++ {
		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}
		z := i / (c_in)
		x := i % (c_in)
		bias[z][x] = complex(value, 0)
	}
	return bias
}

func encoding_filter_bias(c_in int, c_out int, filters int, filter_size int, H_out int, W_out int, filter_value [][][]complex128, bias_value [][]complex128, encoder ckks.Encoder, params ckks.Parameters) (pt_w [][]*ckks.Plaintext, pt_b [][]*ckks.Plaintext) {
	enc_filter := make([][][]complex128, c_out)
	enc_bias := make([][][]complex128, c_out)
	for m := 0; m < c_out; m++ {
		enc_filter[m] = make([][]complex128, c_in)
		enc_bias[m] = make([][]complex128, c_in)
		for l := 0; l < c_in; l++ {
			enc_filter[m][l] = make([]complex128, H_out*W_out*filters)
			enc_bias[m][l] = make([]complex128, H_out*W_out*filters)
			for i := 0; i < H_out; i++ {
				for j := 0; j < W_out; j++ {
					for k := 0; k < filter_size; k++ {
						for n := 0; n < filter_size; n++ {
							enc_filter[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = filter_value[m][l][k*3+n]
							if k == 0 && n == 0 {
								enc_bias[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = bias_value[m][l]
							} else {
								enc_bias[m][l][(W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
							}
						}
					}
				}
			}
		}
	}

	pt_w = make([][]*ckks.Plaintext, c_out)
	pt_b = make([][]*ckks.Plaintext, c_out)
	for m := 0; m < c_out; m++ {
		// fmt.Println(m)
		pt_w[m] = make([]*ckks.Plaintext, c_in)
		pt_b[m] = make([]*ckks.Plaintext, c_in)
		for l := 0; l < c_in; l++ {
			// fmt.Println(l)
			pt_w[m][l] = encoder.EncodeNew(enc_filter[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			pt_b[m][l] = encoder.EncodeNew(enc_bias[m][l], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}
	return pt_w, pt_b
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

	// Scheme context and keys
	//Preparing data
	sample_num := 1
	// x, y := ReadCsv("mnist_test.csv", sample_num)
	x := make([][]complex128, sample_num)
	x[0] = make([]complex128, 784)
	for i := 0; i < len(x[0]); i++ {
		x[0][i] = 1
	}
	fmt.Println(len(x))
	fmt.Println(len(x[0]))
	// fmt.Println(len(y))

	F_H := 28
	F_W := 28
	party1_cols := 3
	//vector<int> party1_index(F_H_out);
	//vector<int> party2_index(F_H*F_W);

	filter_size := 3
	stride := 1
	F_H_out := (F_H - filter_size + 1) / stride
	F_W_out := (F_W - filter_size + 1) / stride
	F_H_out_1 := (F_H_out - filter_size + 1) / stride
	F_W_out_1 := (F_W_out - filter_size + 1) / stride
	filters := filter_size * filter_size

	features1 := make([][]complex128, sample_num)
	features2 := make([][]complex128, sample_num)
	for i := 0; i < sample_num; i++ {
		features1[i] = make([]complex128, F_H_out*F_W_out*filters)
		features2[i] = make([]complex128, F_H_out*F_W_out*filters)
	}

	//Setting filter 0
	c_in_0 := 1
	c_out_0 := 32
	c_in_1 := 32
	c_out_1 := 64
	filters0_value := Readcsv_filter("00.csv", filter_size, filter_size, c_in_0, c_out_0)
	bias0_value := Readcsv_bias("01.csv", c_in_0, c_out_0)
	filters1_value := Readcsv_filter("10.csv", filter_size, filter_size, c_in_1, c_out_1)
	bias1_value := Readcsv_bias("11.csv", c_in_1, c_out_1)

	// pt_w0, pt_b0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, filters0_value, bias0_value, encoder, params)
	// pt_w1, pt_b1 := encoding_filter_bias(c_in_1, c_out_1, filters, filter_size, F_H_out_1, F_W_out_1, filters1_value, bias1_value, encoder, params)
	fmt.Printf("Read filter and bias")
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

	rots := []int{1, 2, 3, 4, 5, 6, 7, 8, F_H_out * F_W_out * filters, 9, 234, 243, 17, 34, 231, 248, 265, 462, 479, 496}
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

	pt_w0, pt_b0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, filters0_value, bias0_value, encoder, params)
	pt_w1, pt_b1 := encoding_filter_bias(c_in_1, c_out_1, filters, filter_size, F_H_out_1, F_W_out_1, filters1_value, bias1_value, encoder, params)
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
	ct_result := make([][]*ckks.Ciphertext, sample_num)
	for i := 0; i < sample_num; i++ {

		evaluator.Add(ct_f1[i], ct_f2[i], ct_f1[i])
		fmt.Println(ct_f1[i].Level())

		t1 := time.Now()
		ct_result[i] = Conv(params, ct_f1, pt_w0, pt_b0, c_in_0, c_out_0, filters, F_H_out*F_W_out*filters, evaluator, decryptor, encoder)
		// printDebug(params, ct_result[i], decryptor, encoder)
		t2 := time.Now()
		fmt.Println("Conv time: ", t2.Sub(t1))
		for j := 0; j < c_out_0; j++ {
			ct_result[i][j] = relu(params, ct_result[i][j], evaluator, decryptor, encoder)
		}
		// fmt.Println(ct_result[i].Level())
		// ct_result[i] = relu(params, ct_result[i], evaluator, decryptor, encoder)
		t3 := time.Now()
		fmt.Println("ReLU time: ", t3.Sub(t2))
		// fmt.Println(ct_result[i].Level())
		for j := 0; j < c_out_0; j++ {
			ct_result[i][j] = avg_pooling(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder)
		}
		t4 := time.Now()
		// fmt.Println(ct_result[i].Level())
		fmt.Println("Pooling time: ", t4.Sub(t3))
		// printDebug(params, ct_result[i], decryptor, encoder)
		for j := 0; j < c_out_0; j++ {
			ct_result[i][j] = rot_for_rep_conv(params, ct_result[i][j], F_H_out*F_W_out*filters, evaluator, decryptor, encoder)
		}
		t5 := time.Now()
		fmt.Println("Replication time: ", t5.Sub(t4))
		// ct_result[i].Rescale()
		// printDebug(params, ct_result[i], decryptor, encoder)
		//conv again
		ct_result[i] = Conv(params, ct_result[i], pt_w1, pt_b1, c_in_0, c_out_0, filters, F_H_out*F_W_out*filters, evaluator, decryptor, encoder)
		// printDebug(params, ct_result[i], decryptor, encoder)

	}

}
func rot_for_rep_conv(params ckks.Parameters, data *ckks.Ciphertext, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (ct *ckks.Ciphertext) {
	rots_for_pool := []int{17, 34, 231, 248, 265, 462, 479, 496}
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

func avg_pooling(params ckks.Parameters, data *ckks.Ciphertext, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (ct *ckks.Ciphertext) {
	rots_for_pool := []int{9, 234, 243}
	dup := evaluator.RotateNew(data, size)
	evaluator.Add(dup, data, dup)
	for i := 0; i < len(rots_for_pool); i++ {
		tmp := evaluator.RotateNew(data, rots_for_pool[i])
		evaluator.Add(dup, tmp, dup)
	}
	mask := make([]complex128, size)
	for i := 0; i < size; i++ {
		if i%18 == 0 {
			mask[i] = complex(1, 0)
		} else {
			mask[i] = complex(0, 0)
		}
	}
	pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	evaluator.Mul(dup, pt_mask, dup)
	evaluator.MultByConst(dup, 0.25, dup)
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
	// levels := make([]int, degree+1)
	// levels[1] = 0
	// levels[0] = 0

	// for i := 2; i <= degree; i++ {
	// 	minlevel := i
	// 	cand := -1
	// 	for j := 1; j <= i/2; j++ {
	// 		k := i - j
	// 		newlevel := Max(levels[j], levels[k]) + 1
	// 		if newlevel < minlevel {
	// 			cand = j
	// 			minlevel = newlevel
	// 		}
	// 	}
	// 	levels[i] = minlevel

	// 	// temp := powers[cand].CopyNew()
	// 	fmt.Println("i", i)
	// 	fmt.Println("cand", cand)
	// 	evaluator.Mul(powers[cand], powers[i-cand], powers[i])
	// 	evaluator.Relinearize(powers[i], powers[i])
	// }

	powers[1] = evaluator.MulNew(powers[0], powers[0])
	evaluator.Relinearize(powers[1], powers[1])

	for i := 1; i <= degree; i++ {
		// fmt.Println(i)
		evaluator.MultByConst(powers[i-1], coeffs[i], powers[i-1])
		// rs[i] = evaluator.MultByConstNew(powers[i], 1)
		// evaluator.Rescale(powers[i], params.DefaultScale(), powers[i])
	}
	for i := 1; i < degree; i++ {
		evaluator.Add(powers[0], powers[i], powers[1])
	}
	ct = evaluator.AddConstNew(powers[1], coeffs[0])
	return ct
}
func relu(params ckks.Parameters, data *ckks.Ciphertext, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result *ckks.Ciphertext) {
	degree := 2
	coeffs := make([]float64, degree+1)
	coeffs[0] = 0.00001
	coeffs[1] = 1.0
	coeffs[2] = 1.0
	result = tree_cipher(params, data, coeffs, degree, evaluator, decryptor, encoder)
	return result
}
func Conv(params ckks.Parameters, data []*ckks.Ciphertext, weights [][]*ckks.Plaintext, bias [][]*ckks.Plaintext, channel_in int, channel_out int, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result_ct []*ckks.Ciphertext) {

	result_ct = make([]*ckks.Ciphertext, channel_out)
	for i := 0; i < channel_out; i++ {
		result_ct[i] = combine_channels(params, data, weights[i], bias[i], channel_in, filters, size, evaluator, decryptor, encoder)
	}
	return result_ct
}

func combine_channels(params ckks.Parameters, data []*ckks.Ciphertext, weights []*ckks.Plaintext, bias []*ckks.Plaintext, channel_in int, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (result *ckks.Ciphertext) {
	result = conv(params, data[0], weights[0], bias[0], filters, size, evaluator, decryptor, encoder)
	for j := 1; j < channel_in; j++ {
		tmp := conv(params, data[j], weights[j], bias[j], filters, size, evaluator, decryptor, encoder)
		evaluator.Add(result, tmp, result)
	}
	return result
}
func conv(params ckks.Parameters, data *ckks.Ciphertext, weights *ckks.Plaintext, bias *ckks.Plaintext, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (dup *ckks.Ciphertext) {
	evaluator.Mul(data, weights, data)
	dup = evaluator.RotateNew(data, size)
	// values := encoder.Decode(decryptor.DecryptNew(dup), params.LogSlots())
	evaluator.Add(dup, data, dup)
	for i := 1; i < filters; i++ {
		tmp := evaluator.RotateNew(data, i)
		evaluator.Add(dup, tmp, dup)
	}
	mask := make([]complex128, size)
	for i := 0; i < size; i++ {
		if i%9 == 0 {
			mask[i] = complex(1, 0)
		} else {
			mask[i] = complex(0, 0)
		}
	}
	pt_mask := encoder.EncodeNew(mask, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	evaluator.Mul(dup, pt_mask, dup)
	evaluator.Add(dup, bias, dup)

	// fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", values[0], values[1], values[2], values[3])
	return dup
}

func printDebug(params ckks.Parameters, ciphertext *ckks.Ciphertext, decryptor ckks.Decryptor, encoder ckks.Encoder) (valuesTest []complex128) {

	valuesTest = encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	// fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	// precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	// fmt.Println(precStats.String())
	fmt.Println(valuesTest[0])
	fmt.Println(valuesTest[18])
	fmt.Println(valuesTest[36])
	fmt.Println(valuesTest[234])

	return
}
