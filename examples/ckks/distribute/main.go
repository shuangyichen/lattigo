package main

import (
	"C"
	"fmt"
	"math/rand"
	"net"
	"os"

	"bufio"
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
var crtgSharesQStr [][][]string
var crtgSharesPStr [][][]string

// func (c *SafeStringArray) Update(str string) {
// 	c.mu.Lock()
// 	c.stringArray[c.numContents] = str
// 	c.numContents++
// 	c.mu.Unlock()
// }
func squeezedArray(input [][]uint64) []uint64 {
	var result = []uint64{}
	for _, arr := range input {
		for _, item := range arr {
			result = append(result, item)
		}
	}
	return result

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

func ServerSetup(serverAddress string, numPeers int) {

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
	// paramSet := bootstrapping.DefaultParametersSparse[2]
	// ckksParams := paramSet.SchemeParams
	// params, err := ckks.NewParametersFromLiteral(ckksParams)
	// ckg := dckks.NewCKGProtocol(params)
	// prng, err := utils.NewKeyedPRNG([]byte{'l', 'a', 't', 't', 'i', 'g', 'o'})
	ckgCombined := ckg.AllocateShare()
	// ringQP, _ := ring.NewRing(params.N(), append(params.Q(), params.P()...))
	pk := ckks.NewPublicKey(params)
	var _ drlwe.CollectivePublicKeyGenerator = dckks.NewCKGProtocol(params)
	// cpkgShares := make([]dckks.CKGShare, numPeers)

	for peerIdx := range clientIPs {
		cpkgShares[peerIdx] = &drlwe.CKGShare{
			Value: ringqp.Poly{},
		}
		// fmtPrintln(clientIPs[peerIdx])
		coeffsQ := polyCoeffsDecode(cpkgSharesQStr[peerIdx])
		coeffsP := polyCoeffsDecode(cpkgSharesPStr[peerIdx])

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
	publicKeyStr = ""
	// pkContent := pk.Value

	//////////////////////////
	// for itemIdx := range pkContent {
	publicKeyStr += polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "&" + polyCoeffsEncode(ckgCombined.Value.P.Coeffs) + "\n"
	// }
	// publicKeyStr += "\n"
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
	galEl := params.GaloisElementForRowRotation()
	rotKeySet := ckks.NewRotationKeySet(params, []uint64{galEl})
	var _ drlwe.RotationKeyGenerator = dckks.NewRotKGProtocol(params)
	for peerIdx := range clientIPs {
		rtgShares[peerIdx] = drlwe.RTGShare{
			Value: [][]ringqp.Poly{},
		}
		polyArr := make([][]ringqp.Poly, len(crtgSharesQStr[peerIdx]))

		for row := range crtgSharesQStr[peerIdx] {
			polyArr[row] = make([]ringqp.Poly, len(crtgSharesQStr[peerIdx][row]))
			for col := range crtgSharesQStr[peerIdx][row] {
				coeffsQ := polyCoeffsDecode(crtgSharesQStr[peerIdx][row][col])
				coeffsP := polyCoeffsDecode(crtgSharesPStr[peerIdx][row][col])

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

}

func ClientSetup(serverAddress string) {

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
	_, err = bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Connected")
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
	p := new(Party)
	p.CKGProtocol = dckks.NewCKGProtocol(params)
	p.s = ckks.NewSecretKey(params) //sk0Shards[0]
	p.s1 = p.AllocateShare()

	crp := p.SampleCRP(prng)

	p.GenShare(p.s, crp, p.s1)
	fmt.Println("Gen pk")

	toSendString := ""
	toSendString += polyCoeffsEncode(p.s1.Value.Q.Coeffs)
	toSendString += "&"
	toSendString += polyCoeffsEncode(p.s1.Value.P.Coeffs)

	toSendString += "\n"
	fmt.Println("Gen sending strings")
	// pk_data, err := p.s1.MarshalBinary()
	// fmt.Println(len(pk_data))
	// n, err := conn.Write(pk_data)
	// fmt.Println(n)
	// if err != nil {
	// 	fmt.Printf("Error writing: %#v\n", err)
	// 	return
	// }
	// n, err := conn.Write([]byte(pk_data))
	// fmt.Println(n)
	conn.Write([]byte(toSendString))
	_, err = bufio.NewReader(conn).ReadString('\n')
	if err != nil {
		panic(err)
	}
	fmt.Println("Receiving CPK")
	////////////////////////Receving CPK//////////////////////////////////////

	//Generating rotation key share
	// p2 := new(PartyRTG)
	// p2.RTGProtocol = dckks.NewRotKGProtocol(params)
	// p2.s = p.s //ckks.NewSecretKey(params)
	// p2.share = p2.AllocateShare()
	// var _ drlwe.RotationKeyGenerator = p2.RTGProtocol
	// crp_rtg := p2.SampleCRP(prng)
	// galEl := params.GaloisElementForRowRotation()
	// p2.GenShare(p2.s, galEl, crp_rtg, p2.share)
	// toSendString = ""
	// for row := range p2.share.Value {
	// 	for col := range p2.share.Value[row] {
	// 		toSendString += polyCoeffsEncode(p2.share.Value[row][col].Q.Coeffs)
	// 		toSendString += "&"
	// 		toSendString += polyCoeffsEncode(p2.share.Value[row][col].P.Coeffs)
	// 		toSendString += "#"
	// 	}
	// 	toSendString += "*"
	// }
	// toSendString += "\n"
	// conn.Write([]byte(toSendString))
	// fmt.Println("Sending rot k share")
	///////////////////////////////////////////////

	//////////////Generating relin key share///////////////////////////
	// p3 := new(PartyRKG)
	// p3.RKGProtocol = dckks.NewRKGProtocol(params)
	// p3.sk = ckks.NewSecretKey(params) // p.s
	// p3.ephSk, p3.share1, p3.share2 = p3.AllocateShare()
	// crp_rkg := p3.SampleCRP(prng)
	// p3.GenShareRoundOne(p3.sk, crp_rkg, p3.ephSk, p3.share1)
	// data, err := p3.share1.MarshalBinary()
	// conn.Write([]byte(data))

	//////////////////////////////////////////////////////////////////

}

func Serverinference() {

}

func Clientinference() {

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
	if appType == "client" {
		inputs := make([]float64, inLen)
		for i := range inputs {
			inputs[i] = rand.Float64()
		}
		fmt.Println("Client setup")
		ClientSetup("localhost:8080")
		// clientPhase2(inputs, cpk, shamirShare, id, "localhost:8080", robust, logDegree, scale, 0.5)
	} else {
		fmt.Println("Server setup")
		ServerSetup("localhost:8080", numPeers)
		// serverPhase2("localhost:8080", numPeers, robust, 0.5, logDegree, scale, inLen)
	}
}
