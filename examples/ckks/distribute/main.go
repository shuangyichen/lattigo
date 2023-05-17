package main

import (
	"C"
	"fmt"
	"math/rand"
	"net"
	"os"

	"bufio"
	"encoding/csv"
	"strconv"

	// "strconv"
	// "strings"
	"sync"

	"github.com/tuneinsight/lattigo/v3/ckks"
	"github.com/tuneinsight/lattigo/v3/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v3/dckks"
	"github.com/tuneinsight/lattigo/v3/drlwe"

	"github.com/tuneinsight/lattigo/v3/ring"
	"github.com/tuneinsight/lattigo/v3/rlwe"
)
import (
	"strings"
	"time"

	//

	"math"
	"math/bits"

	"github.com/tuneinsight/lattigo/v3/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v3/utils"
)

type SafeStringArray struct {
	mu          sync.Mutex
	stringArray []string
	numContents int
}

var publicKeyStr string
var cpkgSharesQStr []string
var cpkgSharesPStr []string
var cpkgShares []*drlwe.CKGShare
var cpkgWrite []sync.Mutex
var pkStrWrite []sync.Mutex
var done []sync.Mutex
var crtgWrite []sync.Mutex
var rtStrWrite []sync.Mutex
var rtdone []sync.Mutex
var rotKeyStr string
var crkgR1Str string
var crtgSharesQStr [][][]string
var crtgSharesPStr [][][]string

var crkgSharesQStr [][][][]string
var crkgSharesPStr [][][][]string
var crkgWrite []sync.Mutex
var crkg2Write []sync.Mutex
var rkgdone []sync.Mutex
var rlkr1StrWrite []sync.Mutex
var ctStr []string
var ctWrite []sync.Mutex
var ct_done []sync.Mutex
var pdWrite []sync.Mutex
var pd_done []sync.Mutex
var pdStr []string

// func (c *SafeStringArray) Update(str string) {
// 	c.mu.Lock()
// 	c.stringArray[c.numContents] = str
// 	c.numContents++
// 	c.mu.Unlock()
// }
func GetMinimumLevelForBootstrapping(lambda int, scale float64, nParties int, moduli []uint64) (minLevel, logBound int, ok bool) {
	logBound = lambda + int(math.Ceil(math.Log2(scale)))
	maxBound := logBound + bits.Len64(uint64(nParties))
	minLevel = -1
	logQ := 0
	for i := 0; logQ < maxBound; i++ {
		if i >= len(moduli) {
			return 0, 0, false
		}

		logQ += bits.Len64(moduli[i])
		minLevel++
	}
	if len(moduli) < minLevel {
		return 0, 0, false
	}

	return minLevel, logBound, true
}

// NewAdditiveShareBigint instantiates a new additive share struct composed of "n" big.Int elements
func NewAdditiveShareBigint(params ckks.Parameters, logSlots int) *rlwe.AdditiveShareBigint {
	dslots := 1 << logSlots
	if params.RingType() == ring.Standard {
		dslots *= 2
	}
	return rlwe.NewAdditiveShareBigint(params.Parameters, dslots)
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

func arrayOfIntToString(a []uint64, delim string) string {
	return strings.Trim(strings.Replace(fmt.Sprint(a), " ", delim, -1), "[]")
	//return strings.Trim(strings.Join(strings.Split(fmt.Sprint(a), " "), delim), "[]")
	//return strings.Trim(strings.Join(strings.Fields(fmt.Sprint(a)), delim), "[]")
}

func squeezedArray(input [][]uint64) []uint64 {
	var result = []uint64{}
	for _, arr := range input {
		for _, item := range arr {
			result = append(result, item)
		}
	}
	return result

}

func unsqueezedArray(input []uint64, num int) [][]uint64 {
	resLoc := make([][]uint64, num)
	length := len(input)
	piece := length / num
	for i := 0; i < num; i++ {
		resLoc[i] = input[piece*i : piece*(i+1)]
	}
	return resLoc
}

func polyCoeffsEncode(coeffs [][]uint64) string {
	res := ""
	for dimZeroCounter := range coeffs {
		// fmt.Println(dimZeroCounter)
		tmp := coeffs[dimZeroCounter]
		tmpRes := arrayToString(tmp)
		res += tmpRes + " "
	}
	return res
}

func arrayToString(arr []uint64) string {
	init := strconv.FormatUint(arr[0], 10)
	// fmt.Println(len(arr))
	for cntr := 1; cntr < len(arr); cntr++ {
		init = init + "." + strconv.FormatUint(arr[cntr], 10)
	}
	return init
}
func polyCoeffsDecode(str string) [][]uint64 {
	//fmt.Println("I am called")
	polyCoeffsArr := strings.Split(str, " ")
	// fmt.Println(len(polyCoeffsArr))
	polyCoeffsArr = polyCoeffsArr[0 : len(polyCoeffsArr)-1]
	//fmt.Println("I am allocating")
	res := make([][]uint64, len(polyCoeffsArr))
	//fmt.Println("I allocated")
	for polyCounter := range res {
		tmp := stringToArray(polyCoeffsArr[polyCounter])
		res[polyCounter] = tmp
	}
	//fmt.Println("I am returning")
	return res
}
func stringToArrayOfUint(strIn string) []uint64 {
	//fmt.Println("string to arr called")
	//fmt.Println(strIn[0:20])
	elements := strings.Split(strIn, ",")
	// fmt.Println("string to arr allocating")
	// fmt.Println(len(elements))
	resLoc := make([]uint64, len(elements))
	//fmt.Println("string to arr allocated")
	for counter := range resLoc {
		resLoc[counter], _ = strconv.ParseUint(elements[counter], 10, 64)
	}
	//fmt.Println("string to arr returned")
	return resLoc
}
func stringToArray(strIn string) []uint64 {
	//fmt.Println("string to arr called")
	//fmt.Println(strIn[0:20])
	elements := strings.Split(strIn, ".")
	// fmt.Println("string to arr allocating")
	// fmt.Println(len(elements))
	resLoc := make([]uint64, len(elements))
	//fmt.Println("string to arr allocated")
	for counter := range resLoc {
		resLoc[counter], _ = strconv.ParseUint(elements[counter], 10, 64)
	}
	//fmt.Println("string to arr returned")
	return resLoc
}

func handleClientrot(connections []net.Conn, idx int) {
	conn := connections[idx]
	message, err := bufio.NewReader(conn).ReadString('\n')
	fmt.Println("Receiving client rot key shares")
	if err != nil {
		panic(err)
	}
	rotsharesArr := strings.Split(message, "*")
	rotsharesArr = rotsharesArr[0 : len(rotsharesArr)-1]
	rows := len(rotsharesArr)
	// fmt.Println(rows) //4
	crtgSharesQStr[idx] = make([][]string, rows)
	crtgSharesPStr[idx] = make([][]string, rows)
	for row := range rotsharesArr {
		rotsharesArrtmp := strings.Split(rotsharesArr[row], "#")
		rotsharesArrtmp = rotsharesArrtmp[0 : len(rotsharesArrtmp)-1]
		crtgSharesQStr[idx][row] = make([]string, len(rotsharesArrtmp))
		crtgSharesPStr[idx][row] = make([]string, len(rotsharesArrtmp))
		// fmt.Println(len(rotsharesArrtmp)) //0
		for col := range rotsharesArrtmp {
			// fmt.Println(col) //0
			rotsharesArrtmp2 := strings.Split(rotsharesArrtmp[col], "&")
			// fmt.Println(len(rotsharesArrtmp2)) //2
			crtgSharesQStr[idx][row][col] = rotsharesArrtmp2[0]
			crtgSharesPStr[idx][row][col] = rotsharesArrtmp2[1]
		}
	}
	crtgWrite[idx].Unlock()

}

func handleClientrrlkr1(connections []net.Conn, idx int) {
	//////
	fmt.Println("Handling client rlk")
	// numPeers := len(connections)
	conn := connections[idx]
	// conn.Write([]byte(strconv.Itoa(idx) + " " + strconv.Itoa(numPeers) + "\n"))
	// fmt.Println("Sending client index")
	///////
	// conn := connections[idx]
	message, err := bufio.NewReader(conn).ReadString('\n')
	fmt.Println("Receiving client relin r1 key shares")
	if err != nil {
		panic(err)
	}
	rlksharesArr := strings.Split(message, "*")
	rlksharesArr = rlksharesArr[0 : len(rlksharesArr)-1]
	rows := len(rlksharesArr)
	// fmt.Println(rows) //4
	crkgSharesQStr[idx] = make([][][]string, rows)
	crkgSharesPStr[idx] = make([][][]string, rows)
	for row := range rlksharesArr {
		rlksharesArrtmp := strings.Split(rlksharesArr[row], "@")
		rlksharesArrtmp = rlksharesArrtmp[0 : len(rlksharesArrtmp)-1]
		crkgSharesQStr[idx][row] = make([][]string, len(rlksharesArrtmp))
		crkgSharesPStr[idx][row] = make([][]string, len(rlksharesArrtmp))
		// fmt.Println(len(rotsharesArrtmp)) //0
		for col := range rlksharesArrtmp {
			rlksharesArrtmp2 := strings.Split(rlksharesArr[row], "#")
			rlksharesArrtmp2 = rlksharesArrtmp2[0 : len(rlksharesArrtmp2)-1]
			crkgSharesQStr[idx][row][col] = make([]string, len(rlksharesArrtmp2))
			crkgSharesPStr[idx][row][col] = make([]string, len(rlksharesArrtmp2))
			for id := range rlksharesArrtmp2 {
				rlksharesArrtmp3 := strings.Split(rlksharesArrtmp2[id], "&")
				crkgSharesQStr[idx][row][col][id] = rlksharesArrtmp3[0]
				crkgSharesPStr[idx][row][col][id] = rlksharesArrtmp3[1]
			}

		}
	}
	crkgWrite[idx].Unlock()
	rlkr1StrWrite[idx].Lock()
	conn.Write([]byte(crkgR1Str))

}

func handleClientrrlkr2(connections []net.Conn, idx int) {

	conn := connections[idx]
	message, err := bufio.NewReader(conn).ReadString('\n')
	fmt.Println("Receiving client relin r2 key shares")
	// fmt.Println("Receiving client rot key shares")
	if err != nil {
		panic(err)
	}
	rlksharesArr := strings.Split(message, "*")
	rlksharesArr = rlksharesArr[0 : len(rlksharesArr)-1]
	rows := len(rlksharesArr)
	// fmt.Println(rows) //4
	crkgSharesQStr[idx] = make([][][]string, rows)
	crkgSharesPStr[idx] = make([][][]string, rows)
	for row := range rlksharesArr {
		rlksharesArrtmp := strings.Split(rlksharesArr[row], "@")
		rlksharesArrtmp = rlksharesArrtmp[0 : len(rlksharesArrtmp)-1]
		crkgSharesQStr[idx][row] = make([][]string, len(rlksharesArrtmp))
		crkgSharesPStr[idx][row] = make([][]string, len(rlksharesArrtmp))
		// fmt.Println(len(rotsharesArrtmp)) //0
		for col := range rlksharesArrtmp {
			rlksharesArrtmp2 := strings.Split(rlksharesArr[row], "#")
			rlksharesArrtmp2 = rlksharesArrtmp2[0 : len(rlksharesArrtmp2)-1]
			crkgSharesQStr[idx][row][col] = make([]string, len(rlksharesArrtmp2))
			crkgSharesPStr[idx][row][col] = make([]string, len(rlksharesArrtmp2))
			for id := range rlksharesArrtmp2 {
				rlksharesArrtmp3 := strings.Split(rlksharesArrtmp2[id], "&")
				crkgSharesQStr[idx][row][col][id] = rlksharesArrtmp3[0]
				crkgSharesPStr[idx][row][col][id] = rlksharesArrtmp3[1]
			}

		}
	}
	crkg2Write[idx].Unlock()
}

func handleClientcpk(connections []net.Conn, idx int) {

	////// transmit the index of the client and number of peers
	fmt.Println("Handling client")
	numPeers := len(connections)
	conn := connections[idx]
	conn.Write([]byte(strconv.Itoa(idx) + " " + strconv.Itoa(numPeers) + "\n"))
	fmt.Println("Sending client index")

	/////// Receive shares and generate collective public key
	message, err := bufio.NewReader(conn).ReadString('\n')

	fmt.Println("Receiving client pk shares")
	if err != nil {
		panic(err)
	}
	sharesArr := strings.Split(message, "&")
	// fmt.Println(len(sharesArr))
	sharesQArr := sharesArr[0]
	sharesPArr := sharesArr[1]

	cpkgSharesQStr[idx] = sharesQArr
	cpkgSharesPStr[idx] = sharesPArr

	// buf := make([]byte, 1)
	// len, err := conn.Read(buf)
	// if err != nil {
	// 	fmt.Printf("Error reading: %#v\n", err)
	// 	return
	// }
	// fmt.Println(len)
	// cpkgShares[idx].UnmarshalBinary(buf[:len])
	// buf[:len]
	cpkgWrite[idx].Unlock()
	pkStrWrite[idx].Lock()
	conn.Write([]byte(publicKeyStr))
	done[idx].Unlock()

}

func ServerSetup(serverAddress string, numPeers int) (pk *rlwe.PublicKey, rotKeySet *rlwe.RotationKeySet, crlk *rlwe.RelinearizationKey) {

	l, err := net.Listen("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error listening:", err.Error())
		os.Exit(1)
	}
	// Close the listener when the application closes.
	defer l.Close()
	// The array to save the address of clients
	clientIPs := make([]net.Conn, numPeers)
	for cntr := 0; cntr < numPeers; cntr++ {
		// Listen for an incoming connection.
		c, err := l.Accept()
		if err != nil {
			fmt.Println("Error connecting:", err.Error())
			return
		}
		// Print client connection address.
		//fmt.Println("Client " + c.RemoteAddr().String() + " connected.")
		// Add connection to list
		clientIPs[cntr] = c
	}
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	params, err := ckks.NewParametersFromLiteral(ckksParams)
	ckg := dckks.NewCKGProtocol(params)
	prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	cpkgShares = make([]*drlwe.CKGShare, numPeers)
	for peerIdx := range cpkgShares {
		cpkgShares[peerIdx] = ckg.AllocateShare()
	}
	cpkgSharesQStr = make([]string, numPeers)
	cpkgSharesPStr = make([]string, numPeers)
	cpkgWrite = make([]sync.Mutex, numPeers)
	done = make([]sync.Mutex, numPeers)
	pkStrWrite = make([]sync.Mutex, numPeers)

	crtgSharesQStr = make([][][]string, numPeers)
	crtgSharesPStr = make([][][]string, numPeers)
	crtgWrite = make([]sync.Mutex, numPeers)
	// rtStrWrite = make([]sync.Mutex, numPeers)
	rtdone = make([]sync.Mutex, numPeers)
	rtgShares := make([]drlwe.RTGShare, numPeers)

	crkgSharesQStr = make([][][][]string, numPeers)
	crkgSharesPStr = make([][][][]string, numPeers)
	crkgWrite = make([]sync.Mutex, numPeers)
	rtStrWrite = make([]sync.Mutex, numPeers)
	rkgdone = make([]sync.Mutex, numPeers)
	crkgWrite = make([]sync.Mutex, numPeers)
	crkg2Write = make([]sync.Mutex, numPeers)
	rkgShares := make([]drlwe.RKGShare, numPeers)
	rlkr1StrWrite = make([]sync.Mutex, numPeers)

	for peerIdx := range cpkgWrite {
		pkStrWrite[peerIdx].Lock()
		cpkgWrite[peerIdx].Lock()
		done[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
		// // rtStrWrite[peerIdx].Lock()
		// rtdone[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientcpk(clientIPs, idx)
	}

	for peerIdx := range cpkgWrite {
		cpkgWrite[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
	}

	//////////////////////////

	ckgCombined := ckg.AllocateShare()
	pk = ckks.NewPublicKey(params)
	var _ drlwe.CollectivePublicKeyGenerator = dckks.NewCKGProtocol(params)
	// cpkgShares := make([]dckks.CKGShare, numPeers)
	q_row := len(ckgCombined.Value.Q.Coeffs)
	p_row := len(ckgCombined.Value.P.Coeffs)
	for peerIdx := range clientIPs {
		cpkgShares[peerIdx] = &drlwe.CKGShare{
			Value: ringqp.Poly{},
		}
		// fmtPrintln(clientIPs[peerIdx])
		coeffsQ := unsqueezedArray(stringToArrayOfUint(cpkgSharesQStr[peerIdx]), q_row) //polyCoeffsDecode(cpkgSharesQStr[peerIdx])
		coeffsP := unsqueezedArray(stringToArrayOfUint(cpkgSharesPStr[peerIdx]), p_row) //polyCoeffsDecode(cpkgSharesPStr[peerIdx])

		poly := ringqp.Poly{
			Q: &ring.Poly{},
			P: &ring.Poly{},
		}
		// poly := ringqp.Poly{P.coeffs == coeffsP, Q.coeffs== coeffsQ}
		poly.P.Coeffs = coeffsP
		poly.Q.Coeffs = coeffsQ
		// ringqp.Poly{{Q:coeffsQ,P:coeffsP}}
		// poly[0] = ring.Poly{Coeffs:coeffsQ}

		cpkgShares[peerIdx].Value = poly

		ckg.AggregateShare(cpkgShares[peerIdx], ckgCombined, ckgCombined)

	}
	crp := ckg.SampleCRP(prng)
	ckg.GenPublicKey(ckgCombined, crp, pk)
	fmt.Println("Gen cpk server side")
	q_array := squeezedArray(ckgCombined.Value.Q.Coeffs)
	p_array := squeezedArray(ckgCombined.Value.P.Coeffs)

	publicKeyStr = ""
	publicKeyStr += arrayOfIntToString(q_array, ",") + "&" + arrayOfIntToString(p_array, ",") + "\n"
	// publicKeyStr += polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "&" + polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "\n"
	fmt.Println("Gen cpk string server side")
	for peerIdx := range cpkgWrite {
		pkStrWrite[peerIdx].Unlock()
	}

	for peerIdx := range done {
		done[peerIdx].Lock()
		// done[peerIdx].Unlock()
	}
	for peerIdx := range crtgWrite {
		crtgWrite[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientrot(clientIPs, idx)
	}
	for peerIdx := range cpkgWrite {
		crtgWrite[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
	}

	rtg := dckks.NewRotKGProtocol(params)
	rtgCombined := rtg.AllocateShare()
	q_row = len(rtgCombined.Value[0][0].Q.Coeffs)
	p_row = len(rtgCombined.Value[0][0].P.Coeffs)
	galEl := params.GaloisElementForRowRotation()
	rotKeySet = ckks.NewRotationKeySet(params, []uint64{galEl})
	var _ drlwe.RotationKeyGenerator = dckks.NewRotKGProtocol(params)
	for peerIdx := range clientIPs {
		rtgShares[peerIdx] = drlwe.RTGShare{
			Value: [][]ringqp.Poly{},
		}
		polyArr := make([][]ringqp.Poly, len(crtgSharesQStr[peerIdx]))

		for row := range crtgSharesQStr[peerIdx] {
			polyArr[row] = make([]ringqp.Poly, len(crtgSharesQStr[peerIdx][row]))
			for col := range crtgSharesQStr[peerIdx][row] {
				coeffsQ := unsqueezedArray(stringToArrayOfUint(crtgSharesQStr[peerIdx][row][col]), q_row)
				coeffsP := unsqueezedArray(stringToArrayOfUint(crtgSharesPStr[peerIdx][row][col]), p_row)

				poly := ringqp.Poly{
					Q: &ring.Poly{},
					P: &ring.Poly{},
				}
				// poly := ringqp.Poly{P.coeffs == coeffsP, Q.coeffs== coeffsQ}
				poly.P.Coeffs = coeffsP
				poly.Q.Coeffs = coeffsQ
				polyArr[row][col] = poly
			}
		}
		rtgShares[peerIdx].Value = polyArr
		// fmt.Println("Rot key share aggregation")
		rtg.AggregateShare(&rtgShares[peerIdx], rtgCombined, rtgCombined)

	}
	crp_rtg := rtg.SampleCRP(prng)
	rtg.GenRotationKey(rtgCombined, crp_rtg, rotKeySet.Keys[galEl])
	fmt.Println("Rot key aggregation")

	// rotKeyStr = ""

	// --------------Relin key Generation -------------------
	for peerIdx := range crkgWrite {
		crkgWrite[peerIdx].Lock()
		rlkr1StrWrite[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientrrlkr1(clientIPs, idx)
	}
	for peerIdx := range crkgWrite {
		crkgWrite[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
	}

	rkg := dckks.NewRKGProtocol(params)
	_, rkg1Combined, rkg2Combined := rkg.AllocateShare()
	q_row_ := len(rkg1Combined.Value[0][0][0].Q.Coeffs)
	p_row_ := len(rkg1Combined.Value[0][0][0].P.Coeffs)
	// galEl := params.GaloisElementForRowRotation()
	// rotKeySet := ckks.NewRotationKeySet(params, []uint64{galEl})
	// var _ drlwe.RotationKeyGenerator = dckks.NewRotKGProtocol(params)
	for peerIdx := range clientIPs {
		rkgShares[peerIdx] = drlwe.RKGShare{
			Value: [][][2]ringqp.Poly{},
		}
		polyArr := make([][][2]ringqp.Poly, len(crkgSharesQStr[peerIdx]))

		for row := range crkgSharesQStr[peerIdx] {
			polyArr[row] = make([][2]ringqp.Poly, len(crkgSharesQStr[peerIdx][row]))
			for col := range crkgSharesQStr[peerIdx][row] {
				// polyArr[row][col] = make([]ringqp.Poly, len(crkgSharesQStr[peerIdx][row][col]))
				for id := range crkgSharesQStr[peerIdx][row][col] {
					coeffsQ := unsqueezedArray(stringToArrayOfUint(crkgSharesQStr[peerIdx][row][col][id]), q_row_)
					coeffsP := unsqueezedArray(stringToArrayOfUint(crkgSharesPStr[peerIdx][row][col][id]), p_row_)
					poly := ringqp.Poly{
						Q: &ring.Poly{},
						P: &ring.Poly{},
					}
					// poly := ringqp.Poly{P.coeffs == coeffsP, Q.coeffs== coeffsQ}
					poly.P.Coeffs = coeffsP
					poly.Q.Coeffs = coeffsQ
					polyArr[row][col][id] = poly
				}
			}
		}
		rkgShares[peerIdx].Value = polyArr
		// fmt.Println("Rot key share aggregation")
		rkg.AggregateShare(&rkgShares[peerIdx], rkg1Combined, rkg1Combined)

	}
	fmt.Println("relin key share r1 aggregation")
	crkgR1Str = ""
	for row := range rkg1Combined.Value {
		for col := range rkg1Combined.Value[row] {
			for idx := range rkg1Combined.Value[row][col] {
				crkgR1Str += arrayOfIntToString(squeezedArray(rkg1Combined.Value[row][col][idx].Q.Coeffs), ",")
				crkgR1Str += "&"
				crkgR1Str += arrayOfIntToString(squeezedArray(rkg1Combined.Value[row][col][idx].P.Coeffs), ",")
				crkgR1Str += "#"
			}
			crkgR1Str += "@"
		}
		crkgR1Str += "*"
	}
	crkgR1Str += "\n"
	for peerIdx := range rlkr1StrWrite {
		rlkr1StrWrite[peerIdx].Unlock()
	}

	for peerIdx := range crkgWrite {
		crkg2Write[peerIdx].Lock()
		// rlkr1StrWrite[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientrrlkr2(clientIPs, idx)
	}
	for peerIdx := range crkgWrite {
		crkg2Write[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
	}
	fmt.Println("r2 aggregation")
	for peerIdx := range clientIPs {
		rkgShares[peerIdx] = drlwe.RKGShare{
			Value: [][][2]ringqp.Poly{},
		}
		polyArr := make([][][2]ringqp.Poly, len(crkgSharesQStr[peerIdx]))

		for row := range crkgSharesQStr[peerIdx] {
			polyArr[row] = make([][2]ringqp.Poly, len(crkgSharesQStr[peerIdx][row]))
			for col := range crkgSharesQStr[peerIdx][row] {
				for id := range crkgSharesQStr[peerIdx][row][col] {
					coeffsQ := unsqueezedArray(stringToArrayOfUint(crkgSharesQStr[peerIdx][row][col][id]), q_row_)
					coeffsP := unsqueezedArray(stringToArrayOfUint(crkgSharesPStr[peerIdx][row][col][id]), p_row_)
					poly := ringqp.Poly{
						Q: &ring.Poly{},
						P: &ring.Poly{},
					}
					poly.P.Coeffs = coeffsP
					poly.Q.Coeffs = coeffsQ
					polyArr[row][col][id] = poly
				}
			}
		}
		rkgShares[peerIdx].Value = polyArr
		// fmt.Println("Rot key share aggregation")
		rkg.AggregateShare(&rkgShares[peerIdx], rkg2Combined, rkg2Combined)

	}
	crlk = ckks.NewRelinearizationKey(params)
	rkg.GenRelinearizationKey(rkg1Combined, rkg2Combined, crlk)
	fmt.Println("Generating relin key")
	toSendString := "1\n"
	for i := range clientIPs {
		conn := clientIPs[i]
		conn.Write([]byte(toSendString))
	}
	// conn.Write([]byte(toSendString))
	// crp_rtg := rtg.SampleCRP(prng)
	// rtg.GenRotationKey(rtgCombined, crp_rtg, rotKeySet.Keys[galEl])
	// fmt.Println("Rot key aggregation")
	return
}

func ClientSetup(serverAddress string) (cpk *rlwe.PublicKey, sk *rlwe.SecretKey, index int) {

	//Build up communication
	// var ringPrime uint64 = 0x10000000001d0001
	// var ringPrimeP uint64 = 0xfffffffffffc001
	///////// Start the client and connect to the server.
	//startTime := time.Now()
	rand.Seed(time.Now().UTC().UnixNano())
	//fmt.Println("Connecting to", "tcp", "server", serverAddress)
	conn, err := net.Dial("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error connecting:", err.Error())
		os.Exit(1)
	}
	defer conn.Close()
	idx_str, err := bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Connected")
	index, err = strconv.Atoi(idx_str)
	// moduli := &ckks.Moduli{Qi: []uint64{ringPrime}, Pi: []uint64{ringPrimeP}}
	// params, err := ckks.NewParametersFromLiteral(logDegree, moduli)
	// params.SetScale(scale)
	// params.SetLogSlots(logDegree - 1)
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	prng, _ := utils.NewKeyedPRNG([]byte{'t', 'e', 's', 't'})
	params, err := ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	type Party struct {
		*dckks.CKGProtocol
		s  *rlwe.SecretKey
		s1 *drlwe.CKGShare
	}
	type PartyRTG struct {
		*dckks.RTGProtocol
		s     *rlwe.SecretKey
		share *drlwe.RTGShare
	}

	type PartyRKG struct {
		*dckks.RKGProtocol
		ephSk  *rlwe.SecretKey
		sk     *rlwe.SecretKey
		share1 *drlwe.RKGShare
		share2 *drlwe.RKGShare
	}
	// minLevel, logBound, ok := dckks.GetMinimumLevelForBootstrapping(128, params.DefaultScale(), numPeers, params.Q())

	//////////////////////////Generating Public key share////////////////////////////////
	pk_start := time.Now()
	p := new(Party)
	p.CKGProtocol = dckks.NewCKGProtocol(params)
	p.s = ckks.NewSecretKey(params) //sk0Shards[0]
	p.s1 = p.AllocateShare()

	crp := p.SampleCRP(prng)

	p.GenShare(p.s, crp, p.s1)
	fmt.Println("Gen pk")

	q_array := squeezedArray(p.s1.Value.Q.Coeffs)
	p_array := squeezedArray(p.s1.Value.P.Coeffs)
	toSendString := arrayOfIntToString(q_array, ",") + "&" + arrayOfIntToString(p_array, ",") + "\n"
	n, err := conn.Write([]byte(toSendString))
	fmt.Println("Sending", n)
	cpk_str, err := bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	cpk_arr := strings.Split(cpk_str, "&")
	q_row := len(p.s1.Value.Q.Coeffs)
	p_row := len(p.s1.Value.P.Coeffs)
	coeffsQ := unsqueezedArray(stringToArrayOfUint(cpk_arr[0]), q_row) //polyCoeffsDecode(cpkgSharesQStr[peerIdx])
	coeffsP := unsqueezedArray(stringToArrayOfUint(cpk_arr[1]), p_row) //polyCoeffsDecode(cpkgSharesPStr[peerIdx])
	poly := ringqp.Poly{
		Q: &ring.Poly{},
		P: &ring.Poly{},
	}
	poly.P.Coeffs = coeffsP
	poly.Q.Coeffs = coeffsQ

	cpk_share := p.AllocateShare()
	cpk_share.Value = poly
	cpk = ckks.NewPublicKey(params)
	p.GenPublicKey(cpk_share, crp, cpk)
	fmt.Println("Receiving target pk")
	dur := time.Since(pk_start).Seconds()
	fmt.Println("CPK gen duration: ", dur)
	// ////////////////////////Receving CPK//////////////////////////////////////
	rtk_start := time.Now()
	// Generating rotation key share
	p2 := new(PartyRTG)
	p2.RTGProtocol = dckks.NewRotKGProtocol(params)
	p2.s = p.s //ckks.NewSecretKey(params)
	p2.share = p2.AllocateShare()
	var _ drlwe.RotationKeyGenerator = p2.RTGProtocol
	crp_rtg := p2.SampleCRP(prng)
	galEl := params.GaloisElementForRowRotation()
	p2.GenShare(p2.s, galEl, crp_rtg, p2.share)
	toSendString = ""
	for row := range p2.share.Value {
		for col := range p2.share.Value[row] {
			toSendString += arrayOfIntToString(squeezedArray(p2.share.Value[row][col].Q.Coeffs), ",")
			toSendString += "&"
			toSendString += arrayOfIntToString(squeezedArray(p2.share.Value[row][col].P.Coeffs), ",")
			toSendString += "#"
		}
		toSendString += "*"
	}
	toSendString += "\n"
	n, err = conn.Write([]byte(toSendString))
	fmt.Println("Sending", n)
	fmt.Println("Sending rot k share")

	dur2 := time.Since(rtk_start).Seconds()
	fmt.Println("CRotK gen duration: ", dur2)
	///////////////////////////////////////////////

	//////////////Generating relin key share///////////////////////////
	rlk_start := time.Now()
	p3 := new(PartyRKG)
	p3.RKGProtocol = dckks.NewRKGProtocol(params)
	p3.sk = p.s //ckks.NewSecretKey(params) // p.s
	p3.ephSk, p3.share1, p3.share2 = p3.AllocateShare()
	crp_rkg := p3.SampleCRP(prng)
	p3.GenShareRoundOne(p3.sk, crp_rkg, p3.ephSk, p3.share1)
	toSendString = ""
	for row := range p3.share1.Value {
		for col := range p3.share1.Value[row] {
			for idx := range p3.share1.Value[row][col] {
				toSendString += arrayOfIntToString(squeezedArray(p3.share1.Value[row][col][idx].Q.Coeffs), ",")
				toSendString += "&"
				toSendString += arrayOfIntToString(squeezedArray(p3.share1.Value[row][col][idx].P.Coeffs), ",")
				toSendString += "#"
			}
			toSendString += "@"
		}
		toSendString += "*"
	}
	toSendString += "\n"
	n, err = conn.Write([]byte(toSendString))
	fmt.Println("Sending", n)
	fmt.Println("Sending relin key r1 share")
	message, err := bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Receiving rlk r1")

	q_row_ := len(p3.share1.Value[0][0][0].Q.Coeffs)
	p_row_ := len(p3.share1.Value[0][0][0].P.Coeffs)
	rlksharesArr := strings.Split(message, "*")
	rlksharesArr = rlksharesArr[0 : len(rlksharesArr)-1]
	rows := len(rlksharesArr)
	// fmt.Println(rows) //4
	crkgr1QStr := make([][][]string, rows)
	crkgr1PStr := make([][][]string, rows)
	polyArr := make([][][2]ringqp.Poly, len(crkgr1QStr))
	for row := range rlksharesArr {
		rlksharesArrtmp := strings.Split(rlksharesArr[row], "@")
		rlksharesArrtmp = rlksharesArrtmp[0 : len(rlksharesArrtmp)-1]
		crkgr1QStr[row] = make([][]string, len(rlksharesArrtmp))
		crkgr1PStr[row] = make([][]string, len(rlksharesArrtmp))
		polyArr[row] = make([][2]ringqp.Poly, len(crkgr1QStr[row]))
		// fmt.Println(len(rotsharesArrtmp)) //0
		for col := range rlksharesArrtmp {
			rlksharesArrtmp2 := strings.Split(rlksharesArr[row], "#")
			rlksharesArrtmp2 = rlksharesArrtmp2[0 : len(rlksharesArrtmp2)-1]
			crkgr1QStr[row][col] = make([]string, len(rlksharesArrtmp2))
			crkgr1PStr[row][col] = make([]string, len(rlksharesArrtmp2))
			for id := range rlksharesArrtmp2 {
				rlksharesArrtmp3 := strings.Split(rlksharesArrtmp2[id], "&")
				crkgr1QStr[row][col][id] = rlksharesArrtmp3[0]
				crkgr1PStr[row][col][id] = rlksharesArrtmp3[1]
				coeffsQ := unsqueezedArray(stringToArrayOfUint(crkgr1QStr[row][col][id]), q_row_)
				coeffsP := unsqueezedArray(stringToArrayOfUint(crkgr1PStr[row][col][id]), p_row_)
				poly := ringqp.Poly{
					Q: &ring.Poly{},
					P: &ring.Poly{},
				}
				// poly := ringqp.Poly{P.coeffs == coeffsP, Q.coeffs== coeffsQ}
				poly.P.Coeffs = coeffsP
				poly.Q.Coeffs = coeffsQ
				polyArr[row][col][id] = poly
			}

		}
	}
	p3.share1.Value = polyArr
	p3.GenShareRoundTwo(p3.ephSk, p3.sk, p3.share1, p3.share2)
	fmt.Println("Generating relin key r2 share")
	toSendString = ""
	for row := range p3.share2.Value {
		for col := range p3.share2.Value[row] {
			for idx := range p3.share2.Value[row][col] {
				toSendString += arrayOfIntToString(squeezedArray(p3.share2.Value[row][col][idx].Q.Coeffs), ",")
				toSendString += "&"
				toSendString += arrayOfIntToString(squeezedArray(p3.share2.Value[row][col][idx].P.Coeffs), ",")
				toSendString += "#"
			}
			toSendString += "@"
		}
		toSendString += "*"
	}
	toSendString += "\n"
	// fmt.Println("Sending relin key r2 share")
	conn.Write([]byte(toSendString))
	cost := time.Since(rlk_start).Seconds()
	fmt.Println("Duration:", cost)
	fmt.Println("Sending relin key r2 share")
	sk = p.s.CopyNew()
	message, err = bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	// data, err := p3.share1.MarshalBinary()
	// conn.Write([]byte(data))
	return
	//////////////////////////////////////////////////////////////////

}

func handleClientct(connections []net.Conn, idx int) {

	conn := connections[idx]

	/////// Receive shares and generate collective public key
	message, err := bufio.NewReader(conn).ReadString('\n')

	fmt.Println("Receiving client ct")
	if err != nil {
		panic(err)
	}

	ctStr[idx] = message
	ctWrite[idx].Unlock()
	// pkStrWrite[idx].Lock()
	// conn.Write([]byte(publicKeyStr))
	ct_done[idx].Unlock()

}

func handleClientpd(connections []net.Conn, idx int) {

	conn := connections[idx]

	/////// Receive shares and generate collective public key
	message, err := bufio.NewReader(conn).ReadString('\n')

	fmt.Println("Receiving client pd")
	if err != nil {
		panic(err)
	}

	pdStr[idx] = message
	pdWrite[idx].Unlock()
	// pkStrWrite[idx].Lock()
	// conn.Write([]byte(publicKeyStr))
	pd_done[idx].Unlock()

}

func Serverinference(cpk *rlwe.PublicKey, cRotk *rlwe.RotationKeySet, cRlk *rlwe.RelinearizationKey, numPeers int, serverAddress string) {
	l, err := net.Listen("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error listening:", err.Error())
		os.Exit(1)
	}
	// Close the listener when the application closes.
	defer l.Close()
	// The array to save the address of clients
	conn_set := make([]net.Conn, numPeers)
	for cntr := 0; cntr < numPeers; cntr++ {
		// Listen for an incoming connection.
		c, err := l.Accept()
		if err != nil {
			fmt.Println("Error connecting:", err.Error())
			return
		}
		// Print client connection address.
		//fmt.Println("Client " + c.RemoteAddr().String() + " connected.")
		// Add connection to list
		conn_set[cntr] = c
	}
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	params, _ := ckks.NewParametersFromLiteral(ckksParams)
	var minLevel int
	var ok bool
	if minLevel, _, ok = GetMinimumLevelForBootstrapping(128, params.DefaultScale(), 2, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
		// t.Skip("Not enough levels to ensure correcness and 128 security")
	}

	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: cRlk, Rtks: cRotk})
	// encoder := ckks.NewEncoder(params)
	//Input info
	// F_H := 28
	// F_W := 28
	// filter_size := 3
	// c_in_0 := 1
	// c_out_0 := 4
	// c_in_1 := 4
	// c_out_1 := 4
	// features_num := 100
	// F_H_out_1 := 11
	// F_W_out_1 := 11
	// filters := filter_size * filter_size
	// stride := 1
	// F_H_out := (F_H - filter_size + 1) / stride
	// F_W_out := (F_W - filter_size + 1) / stride
	// rots := []int{1, 2, 3, 4, 5, 6, 7, 8, F_H_out * F_W_out * filters, 9, 234, 243, 17, 34, 465, 482, 499, 930, 947, 964, 198, 216, 36, 234, 252, 18, 468, 486, 936}
	// filters0_value := Readcsv_filter("00.csv", filter_size, filter_size, c_in_0, c_out_0)
	// fmt.Println(filters0_value[0][0])
	// fmt.Println(filters0_value[1][0])
	// fmt.Println(filters0_value[2][0])
	// fmt.Println(filters0_value[3][0])

	// bias0_value := Readcsv_bias("01.csv", c_in_0, c_out_0)
	// filters1_value := Readcsv_filter("30.csv", filter_size, filter_size, c_in_1, c_out_1)
	// bias1_value := Readcsv_bias("31.csv", c_in_1, c_out_1)
	// fc_value := Readcsv_fc("70.csv", 10, features_num)

	// fc_bias := Readcsv_fc_bias("71.csv", 10)

	// pt_w0 := encoding_filter_bias(c_in_0, c_out_0, filters, filter_size, F_H_out, F_W_out, F_H_out+2, F_W_out+2, F_H_out*F_W_out*filters, filters0_value, encoder, params)

	// fmt.Printf("Encoding finished")
	// filters_1 := 18
	// pt_w1 := encoding_filter_bias(c_in_1, c_out_1, filters_1, filter_size, F_H_out_1, F_W_out_1, F_H_out_1+2, F_W_out_1+2, F_H_out*F_W_out*filters, filters1_value, encoder, params)
	// pt_fc := encoding_fc(10, 4, 5, 5, 26*26*9, 100, fc_value, encoder, params)
	fmt.Printf("Encoding finished \n")
	ctStr = make([]string, numPeers)
	ctWrite = make([]sync.Mutex, numPeers)
	ct_done = make([]sync.Mutex, numPeers)
	pdWrite = make([]sync.Mutex, numPeers)
	pd_done = make([]sync.Mutex, numPeers)
	pdStr = make([]string, numPeers)
	for peerIdx := range crkgWrite {
		ctWrite[peerIdx].Lock()
		ct_done[peerIdx].Lock()
		// rlkr1StrWrite[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientct(conn_set, idx)
	}
	for peerIdx := range crkgWrite {
		ctWrite[peerIdx].Lock()
		// crtgWrite[peerIdx].Lock()
	}
	//Recovering ciphertexts
	cttmp := ckks.NewCiphertext(params, params.LogSlots(), params.MaxLevel(), params.DefaultScale())
	col := len(cttmp.Ciphertext.Value[0].Coeffs)
	ct_Arr := make([]*ckks.Ciphertext, numPeers)

	for peerIdx := range ctStr {
		ct_Arr[peerIdx] = ckks.NewCiphertext(params, params.LogSlots(), params.MaxLevel(), params.DefaultScale())
		polyArr := strings.Split(ctStr[peerIdx], "&")
		polyArr = polyArr[:len(polyArr)-1]
		// num_poly := len(polyArr)

		for idx_poly := range polyArr {
			arr := unsqueezedArray(stringToArrayOfUint(polyArr[idx_poly]), col)
			ct_Arr[peerIdx].Value[idx_poly].Coeffs = arr
		}
	}

	//Agggregation
	evaluator.Add(ct_Arr[0], ct_Arr[1], ct_Arr[0])
	fmt.Println("Form ct of complete dataset")

	fmt.Println("Distribute decryption")

	ct_str := ""
	for i := range ct_Arr[0].Ciphertext.Value {
		ct_str += arrayOfIntToString(squeezedArray(ct_Arr[0].Ciphertext.Value[i].Coeffs), ",")
		ct_str += "&"
	}
	ct_str += "\n"
	for idx := 0; idx < numPeers; idx++ {
		conn := conn_set[idx]
		_, err := conn.Write([]byte(ct_str))
		if err != nil {
			panic(err)
		}
	}

	for peerIdx := range crkgWrite {
		pdWrite[peerIdx].Lock()
		pd_done[peerIdx].Lock()
	}
	for idx := 0; idx < numPeers; idx++ {
		go handleClientpd(conn_set, idx)
	}
	for peerIdx := range crkgWrite {
		pdWrite[peerIdx].Lock()
	}
	cks := dckks.NewE2SProtocol(params, 3.19)
	pdtmp := cks.AllocateShare(minLevel)
	col = len(pdtmp.Value.Coeffs)
	pd_Arr := make([]*drlwe.CKSShare, numPeers)

	for peerIdx := range ctStr {
		pd_Arr[peerIdx] = cks.AllocateShare(minLevel)
		arr := unsqueezedArray(stringToArrayOfUint(pdStr[peerIdx]), col)
		pd_Arr[peerIdx].Value.Coeffs = arr
		cks.AggregateShare(pd_Arr[peerIdx], pdtmp, pdtmp)
	}
	secretShare := NewAdditiveShareBigint(params, params.LogSlots())
	cks.GetShare(secretShare, pdtmp, params.LogSlots(), ct_Arr[0], secretShare)
	rec := NewAdditiveShareBigint(params, params.LogSlots())
	// 	for _, p := range P {
	a := rec.Value
	b := secretShare.Value

	for i := range a {
		a[i].Add(a[i], b[i])
	}
	pt := ckks.NewPlaintext(params, ct_Arr[0].Level(), ct_Arr[0].Scale)
	pt.Value.IsNTT = false
	encoder := ckks.NewEncoder(params)
	params.RingQ().SetCoefficientsBigintLvl(pt.Level(), rec.Value, pt.Value)
	_ = encoder.DecodePublic(pt, params.LogSlots(), 0)
	fmt.Println("Decrypt!")
	// 	}
	// ct_final := ckks.NewCiphertext(params, params.LogSlots(), params.MaxLevel(), params.DefaultScale())
	// cks.KeySwitch(ct_Arr[0], pdtmp, ct_final)
	// decryptor.Decrypt(ct_final, ptres)

}

func Clientinference(cpk *rlwe.PublicKey, sk *rlwe.SecretKey, idx int, serverAddress string) {

	rand.Seed(time.Now().UTC().UnixNano())
	//fmt.Println("Connecting to", "tcp", "server", serverAddress)
	conn, err := net.Dial("tcp", serverAddress)
	if err != nil {
		fmt.Println("Error connecting:", err.Error())
		os.Exit(1)
	}
	defer conn.Close()
	paramSet := bootstrapping.DefaultParametersSparse[2]
	ckksParams := paramSet.SchemeParams
	params, _ := ckks.NewParametersFromLiteral(ckksParams)
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, cpk)
	var minLevel, logBound int
	var ok bool
	if minLevel, logBound, ok = GetMinimumLevelForBootstrapping(128, params.DefaultScale(), 2, params.Q()); ok != true || minLevel+1 > params.MaxLevel() {
		// t.Skip("Not enough levels to ensure correcness and 128 security")
	}
	// kgen := ckks.NewKeyGenerator(params)
	///Info of first layer
	F_H := 28
	F_W := 28
	filter_size := 3
	filters := filter_size * filter_size
	stride := 1
	F_H_out := (F_H - filter_size + 1) / stride
	F_W_out := (F_W - filter_size + 1) / stride

	sample_num := 1
	xs, _ := ReadCsv("data/mnist_test.csv", sample_num)
	x := make([][]complex128, sample_num)
	x[0] = xs[0]
	// label := ys[0]
	party1_cols := 14
	features := make([][]complex128, sample_num)
	for i := 0; i < sample_num; i++ {
		features[i] = make([]complex128, F_H_out*F_W_out*filters)
	}

	for i := 0; i < F_H_out; i++ {
		for j := 0; j < F_W_out; j++ {
			for k := 0; k < filter_size; k++ {
				for n := 0; n < filter_size; n++ {
					if idx == 0 {
						if j+n < party1_cols {
							features[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = x[0][(i+k)*F_H+j+n]
							// features2[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
						} else {
							// features2[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = x[0][(i+k)*F_H+j+n]
							features[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
						}
					} else {
						if j+n > party1_cols {
							features[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = x[0][(i+k)*F_H+j+n]
						} else {
							features[0][(F_W_out*(i)+j)*filters+k*filter_size+n] = complex(0, 0)
						}
					}
				}
			}
		}

	}

	pt_f1 := make([]*ckks.Plaintext, sample_num)
	ct_f1 := make([]*ckks.Ciphertext, sample_num)
	for i := 0; i < sample_num; i++ {
		pt_f1[i] = encoder.EncodeNew(features[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ct_f1[i] = encryptor.EncryptNew(pt_f1[i])

	}

	// Ciphertext to string
	ct_str := ""
	for i := range ct_f1[0].Ciphertext.Value {
		ct_str += arrayOfIntToString(squeezedArray(ct_f1[0].Ciphertext.Value[i].Coeffs), ",")
		ct_str += "&"
	}
	ct_str += "\n"
	n, err := conn.Write([]byte(ct_str))
	if err != nil {
		panic(err)
	}

	fmt.Println("Send ct of sub features", n)

	message, err := bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Receiving the ciphertext to decrypt")

	dec_ct := ckks.NewCiphertext(params, params.LogSlots(), params.MaxLevel(), params.DefaultScale())
	col := len(dec_ct.Ciphertext.Value[0].Coeffs)
	polyArr := strings.Split(message, "&")
	polyArr = polyArr[:len(polyArr)-1]
	for idx_poly := range polyArr {
		arr := unsqueezedArray(stringToArrayOfUint(polyArr[idx_poly]), col)
		dec_ct.Value[idx_poly].Coeffs = arr
	}
	cks := dckks.NewE2SProtocol(params, 3.19)

	pd_share := cks.AllocateShare(minLevel)
	secretShare := NewAdditiveShareBigint(params, params.LogSlots())
	cks.GenShare(sk, logBound, params.LogSlots(), dec_ct.Ciphertext.Value[1], secretShare, pd_share)
	pd_array := squeezedArray(pd_share.Value.Coeffs)
	pd_str := arrayOfIntToString(pd_array, ",") + "\n"
	n, err = conn.Write([]byte(pd_str))
	// fmt.Println("Sending", n)
}

func main() {
	args := os.Args[1:]
	appType := args[0]
	fmt.Println(appType)
	// scale := float64(1 << 40)
	// robust := false
	// logDegree := uint64(13)
	numPeers := 2
	inLen := 20000
	server_addr := "10.30.8.11:5000" //"localhost:8080"
	if appType == "client" {
		inputs := make([]float64, inLen)
		for i := range inputs {
			inputs[i] = rand.Float64()
		}
		fmt.Println("Client setup")
		// cpk, sk, idx := ClientSetup(server_addr)
		_, _, _ = ClientSetup(server_addr)
		// time.Sleep(3 * time.Second)
		// Clientinference(cpk, sk, idx, server_addr)
		// clientPhase2(inputs, cpk, shamirShare, id, "localhost:8080", robust, logDegree, scale, 0.5)
	} else {
		fmt.Println("Server setup")
		// cpk, cRotk, cRlk := ServerSetup(server_addr, numPeers)
		_, _, _ = ServerSetup(server_addr, numPeers)
		// Serverinference(cpk, cRotk, cRlk, numPeers, server_addr)
		// serverPhase2("localhost:8080", numPeers, robust, 0.5, logDegree, scale, inLen)
	}
}
