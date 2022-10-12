package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"

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
				y[i] = value //complex(float64(value), 0)
			} else {
				value, err := strconv.Atoi(Rec[i][j])
				// _, err := strconv.Atoi(Rec[i][j])
				if err != nil {
					fmt.Println(err)
				}
				x[i][j-1] = complex(float64(value)/float64(255), 0)

			}
		}
	}
	// 	//Rec, err := ReadCsv.Read()
	// 	value_x0, err := strconv.ParseFloat(Rec[i][0], 64)
	// 	if err != nil {
	// 		fmt.Println(err)
	// 	}
	// 	value_x1, err := strconv.ParseFloat(Rec[i][1], 64)
	// 	if err != nil {
	// 		fmt.Println(err)
	// 	}
	// 	value_x2, err := strconv.ParseFloat(Rec[i][2], 64)
	// 	if err != nil {
	// 		fmt.Println(err)
	// 	}
	// 	value_y, err := strconv.ParseFloat(Rec[i][3], 64)
	// 	if err != nil {
	// 		fmt.Println(err)
	// 	}
	// 	x[i] = make([]complex128)
	// 	x[i][0] = complex(value_x0, 0)
	// 	x[i][1] = complex(value_x1, 0)
	// 	x[i][2] = complex(value_x2, 0)
	// 	y[i] = complex(value_y, 0)
	// }
	return x, y
}

func Read_weights(filename string, output int, fea_num int) (x [][]complex128) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	x = make([][]complex128, output)
	for i := 0; i < output; i++ {
		x[i] = make([]complex128, fea_num)
	}
	// fmt.Println("Sample number", num)

	for i := 0; i < len(Rec); i++ {
		// /x[i] = make([]complex128, fea_num)
		// for j := 0; j < 784; j++ {

		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}
		k := i / 784
		l := i % 784
		x[k][l] = complex(float64(value), 0)

		// }
		// }
	}
	return x
}
func Read_bias(filename string, output int) (x []float64) {
	opencast, err := os.Open(filename)
	if err != nil {
		fmt.Println(err)
	}
	defer opencast.Close()

	ReadCsv := csv.NewReader(opencast)
	Rec, err := ReadCsv.ReadAll()
	x = make([]float64, output)
	// for i := 0; i < output; i++ {
	// 	x[i] = make([]complex128, fea_num)
	// }
	// fmt.Println("Sample number", num)

	for i := 0; i < len(Rec); i++ {
		// /x[i] = make([]complex128, fea_num)
		// for j := 0; j < 784; j++ {

		value, err := strconv.ParseFloat(Rec[i][0], 64)
		if err != nil {
			fmt.Println(err)
		}
		// k := i / 784
		// l := i % 784
		x[i] = float64(value) //complex(float64(value), 0)

		// }
		// }
	}
	return x
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

	rots := []int{1, 784}
	// rots_
	// rots := make([]int, filters+1)
	// for i := 0; i < filters; i++ {
	// 	rots[i] = i
	// }
	// rots[filters] = F_H_out * F_W_out * filters

	kgen = ckks.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPair()
	rlk := kgen.GenRelinearizationKey(sk, 2)
	rotk := kgen.GenRotationKeysForRotations(append(rots, 0-784), true, sk)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotk})
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	encoder = ckks.NewEncoder(params)
	decryptor = ckks.NewDecryptor(params, sk)
	encryptor = ckks.NewEncryptor(params, pk)

	weights_value := Read_weights("10.csv", 10, 784)
	bias_value := Read_bias("11.csv", 10)

	pt_w := make([]*ckks.Plaintext, 10)
	for i := 0; i < 10; i++ {
		pt_w[i] = encoder.EncodeNew(weights_value[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	}

	prediction_correct := 0
	sample_num := 1
	samples := 50
	xs, y := ReadCsv("mnist_test.csv", samples)
	for m := 0; m < 10; m++ {
		fmt.Println("Processing ", m, "-th sample")
		x := make([][]complex128, 1)
		x[0] = xs[m]
		label := y[m]

		// x := make([][]complex128, sample_num)
		// x[0] = make([]complex128, 784)
		// for i := 0; i < len(x[0]); i++ {
		// 	x[0][i] = 1
		// }
		// fmt.Println(len(x))
		// fmt.Println(len(x[0]))
		// fmt.Println(len(y))

		total_feature := 784
		party1_feature_num := 3
		//vector<int> party1_index(F_H_out);
		//vector<int> party2_index(F_H*F_W);

		features1 := make([][]complex128, sample_num)
		features2 := make([][]complex128, sample_num)
		for i := 0; i < sample_num; i++ {
			features1[i] = make([]complex128, total_feature)
			features2[i] = make([]complex128, total_feature)
		}
		// weights_value := make([]complex128, total_feature)

		for i := 0; i < total_feature; i++ {
			if i < party1_feature_num {
				features1[0][i] = x[0][i]
				features2[0][i] = complex(0, 0)
			} else {
				features2[0][i] = x[0][i]
				features1[0][i] = complex(0, 0)
			}
		}
		// weights_value := Read_weights("10.csv", 10, 784)
		// bias_value := Read_bias("11.csv", 10)

		//generating filter
		// for i := 0; i < total_feature; i++ {
		// 	// filters_value[i] = complex(rand.Float64(), 0)
		// 	weights_value[i] = complex(1.0, 0)
		// }

		//encoding filter
		// enc_filter := make([]complex128, F_H_out*F_W_out*filters)
		// for i := 0; i < F_H_out; i++ {
		// 	for j := 0; j < F_W_out; j++ {
		// 		for k := 0; k < filter_size; k++ {
		// 			for n := 0; n < filter_size; n++ {
		// 				enc_filter[(F_W_out*(i)+j)*filters+k*filter_size+n] = filters_value[k*3+n]
		// 			}
		// 		}
		// 	}
		// }
		// rots := []int{1, 784}
		// rots_
		// rots := make([]int, filters+1)
		// for i := 0; i < filters; i++ {
		// 	rots[i] = i
		// }
		// rots[filters] = F_H_out * F_W_out * filters

		// kgen = ckks.NewKeyGenerator(params)

		// sk, pk = kgen.GenKeyPair()
		// rlk := kgen.GenRelinearizationKey(sk, 2)
		// rotk := kgen.GenRotationKeysForRotations(append(rots, 0-784), true, sk)
		// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rotk})
		// // evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

		// encoder = ckks.NewEncoder(params)
		// decryptor = ckks.NewDecryptor(params, sk)
		// encryptor = ckks.NewEncryptor(params, pk)
		// pt_w := make([]*ckks.Plaintext, 10)
		// for i := 0; i < 10; i++ {
		// 	pt_w[i] = encoder.EncodeNew(weights_value[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		// }
		// pt_b := encoder.EncodeNew(bias_value, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

		//encode features
		pt_f1 := make([]*ckks.Plaintext, sample_num)
		pt_f2 := make([]*ckks.Plaintext, sample_num)
		ct_f1 := make([]*ckks.Ciphertext, sample_num)
		ct_f2 := make([]*ckks.Ciphertext, sample_num)
		for i := 0; i < sample_num; i++ {
			// pt_f1[i] = ckks.NewPlaintext(params, params.LogSlots(), params.DefaultScale())
			// pt_f2[i] = ckks.NewPlaintext(params, params.LogSlots(), params.DefaultScale())
			// ct_f1[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
			// ct_f2[i] = ckks.NewCiphertext(params, 1, params.MaxLevel(), params.DefaultScale())
			pt_f1[i] = encoder.EncodeNew(features1[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			pt_f2[i] = encoder.EncodeNew(features2[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			ct_f1[i] = encryptor.EncryptNew(pt_f1[i])
			ct_f2[i] = encryptor.EncryptNew(pt_f2[i])
		}
		//aggregating clients features
		// ct_result := make([]*ckks.Ciphertext, sample_num)
		for i := 0; i < sample_num; i++ {
			evaluator.Add(ct_f1[i], ct_f2[i], ct_f1[i])
			// printDebug(params, ct_f1[i], decryptor, encoder)
			// fmt.Println(ct_f1[i].Level())
			res_ct := FC(params, ct_f1[i], pt_w, bias_value, 784, 784, evaluator, decryptor, encoder)
			prediction_res := Decrypt(params, res_ct, decryptor, encoder)
			fmt.Println("Prediction result: ", prediction_res)
			fmt.Println("Label: ", label)
			if prediction_res == label {
				prediction_correct += 1
				// fmt.Println("Prediction result: ", prediction_res)
				// fmt.Println("Label: ", label)
				fmt.Println("Correct Prediction:", prediction_correct)
			}
			// printDebug(params, ct_result[i], decryptor, encoder)

		}
	}

}

func FC(params ckks.Parameters, data *ckks.Ciphertext, weights []*ckks.Plaintext, bias_value []float64, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (res_ct []*ckks.Ciphertext) {
	res_ct = make([]*ckks.Ciphertext, 10)
	for i := 0; i < 10; i++ {
		fmt.Println("Computing ", i, "-th channel")
		res_ct[i] = fc(params, data, weights[i], filters, size, evaluator, decryptor, encoder)
		evaluator.AddConst(res_ct[i], bias_value[i], res_ct[i])
		// printDebug(params, res_ct[i], decryptor, encoder)
	}
	return res_ct
	// printDebug(params,res_ct[i],decryptor,encoder)
}
func fc(params ckks.Parameters, data *ckks.Ciphertext, weights *ckks.Plaintext, filters int, size int, evaluator ckks.Evaluator, decryptor ckks.Decryptor, encoder ckks.Encoder) (dup *ckks.Ciphertext) {
	res := evaluator.MulNew(data, weights)
	evaluator.Rescale(res, params.DefaultScale(), res)
	dup = evaluator.RotateNew(res, 784)
	// printDebug(params, dup, decryptor, encoder)
	// values := encoder.Decode(decryptor.DecryptNew(dup), params.LogSlots())
	evaluator.Add(dup, data, dup)
	for i := 1; i < 784; i++ {
		evaluator.Rotate(res, 1, res)
		evaluator.Add(dup, res, dup)
	}
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

	// fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", values[0], values[1], values[2], values[3])
	return dup
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
	// fmt.Println(valuesTest[18])
	// fmt.Println(valuesTest[36])
	// fmt.Println(valuesTest[234])

	return
}
