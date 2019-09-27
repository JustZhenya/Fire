#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <bitset>

using namespace std;

// int - степень, int - множитель
typedef bitset<32> Poly;

const size_t buf_size = 100 * 1024 * 1024;
map<int, vector<Poly>> gxes = {
	{2, {0b111}},
	{3, {0b1011}},
	{4, {0b10011, 0b11111}},
	{5, {0b100101}},
	{6, {0b1000011, 0b1001001, 0b1101111}},
	{7, {0b11100111, 0b10000011, 0b10011101}}
};

template<std::size_t N>
void reverse(std::bitset<N>& b) {
	for (std::size_t i = 0; i < N / 2; ++i) {
		bool t = b[i];
		b[i] = b[N - i - 1];
		b[N - i - 1] = t;
	}
}

int NOD(int n1, int n2)
{
	int div;
	if (n1 == n2)  return n1;
	int d = n1 - n2;
	if (d < 0) {
		d = -d;
		div = NOD(n1, d);
	}
	else
		div = NOD(n2, d);
	return div;
}
// Наименьшее общее кратное
int NOK(int n1, int n2)
{
	return n1 * n2 / NOD(n1, n2);
}

bool isPrime(int n)
{
	// Corner case 
	if (n <= 1)
		return false;

	// Check from 2 to n-1 
	for (int i = 2; i < n; i++)
		if (n % i == 0)
			return false;

	return true;
}

class MyPoly
{
public:
	Poly p;

	MyPoly() : p()
	{

	}

	MyPoly(const Poly& _p): p(_p)
	{

	}

	MyPoly operator*(const MyPoly& other)
	{
		MyPoly result;
		for (size_t i = 0; i < other.p.size(); ++i)
		{
			for (size_t j = 0; j < p.size(); ++j)
			{
				if(other.p[i] && p[j])
					result.p.set(i + j);
			}
		}
		return result;
	}

	/*
	def divide(val1,val2):
	a = toList(val1)
	b = toList(val2)
	working=toString(val1)+"\n"

	res=""
	addspace=""

	while len(b) <= len(a) and a:
    		if a[0] == 1:
        		del a[0]
        		for j in range(len(b)-1):
            			a[j] ^= b[j+1]
        		if (len(a)>0):
				working +=addspace+toString(b)+"\n"
				working +=addspace+"-" * (len(b))+"\n"
				addspace+=" "
				working +=addspace+toString(a)+"\n"
				res+= "1"

		else:
        		del a[0]
			working +=addspace+"0" * (len(b))+"\n"
			working +=addspace+"-" * (len(b))+"\n"
			addspace+=" "
			working +=addspace+toString(a)+"\n"

        		res+="0"


	print "Result is\t",res
	print "Remainder is\t",toString(a)

	print "Working is\t\n\n",res.rjust(len(val1)),"\n",
	print "-" * (len(val1)),"\n",working

	return toString(a)
	*/

	MyPoly operator+(const MyPoly& other)
	{
		return other.p | p;
	}

	MyPoly operator/(const MyPoly& divisor)
	{
		uint32_t a_size = get_msbit(), b_size = divisor.get_msbit();

		vector<bool> a, b, result;
		for (size_t i = 0; i < a_size; ++i)
			a.push_back(p[i]);
		for (size_t i = 0; i < b_size; ++i)
			b.push_back(divisor.p[i]);

		reverse(a.begin(), a.end());
		reverse(b.begin(), b.end());

		size_t append_zeroes = b.size() - 1;
		for (size_t i = 0; i < append_zeroes; ++i)
			a.push_back(false);

		while (b.size() <= a.size() && !a.empty())
		{
			if (a[0])
			{
				a.erase(a.begin());
				for (size_t j = 0; j < b.size() - 1; ++j)
					a[j] = (a[j] != b[j + 1]); // xor

				if (!a.empty())
					result.push_back(true);
			}
			else
			{
				a.erase(a.begin());
				result.push_back(false);
			}
		}

		reverse(result.begin(), result.end());

		Poly p;
		for (size_t i = append_zeroes; i < result.size(); ++i)
			p[i - append_zeroes] = result[i];

		return p;
	}

	MyPoly operator%(const MyPoly& divisor)
	{
		uint32_t a_size = get_msbit(), b_size = divisor.get_msbit();

		vector<bool> a, b, result;
		for (size_t i = 0; i < a_size; ++i)
			a.push_back(p[i]);
		for (size_t i = 0; i < b_size; ++i)
			b.push_back(divisor.p[i]);

		reverse(a.begin(), a.end());
		reverse(b.begin(), b.end());

		size_t append_zeroes = b.size() - 1;
		for (size_t i = 0; i < append_zeroes; ++i)
			a.push_back(false);

		while (b.size() <= a.size() && !a.empty())
		{
			if (a[0])
			{
				a.erase(a.begin());
				for (size_t j = 0; j < b.size() - 1; ++j)
					a[j] = (a[j] != b[j + 1]); // xor

				if (!a.empty())
					result.push_back(true);
			}
			else
			{
				a.erase(a.begin());
				result.push_back(false);
			}
		}

		reverse(a.begin(), a.end());

		Poly p;
		for (size_t i = 0; i < a.size(); ++i)
			p[i] = a[i];

		return p;
	}

	/*MyPoly operator%(const MyPoly& divisor)
	{
		size_t a_size = get_msbit(), b_size = divisor.get_msbit();

		size_t result_idx = 0;
		Poly a = p, b = divisor.p, result;

		size_t append_zeroes = b_size - 1;
		a <<= append_zeroes;
		a_size += append_zeroes;

		reverse(a);
		reverse(b);

		while (b_size <= a_size && a_size)
		{
			if (a[a_size - 1])
			{
				//a >>= 1;
				--a_size;

				for (size_t j = 0; j > b_size - 1; --j)
					a[j] = (a[j] != b[j + 1]); // xor

				if (a_size)
					result[result_idx++] = true;
			}
			else
			{
				//a >>= 1;
				--a_size;

				if (a_size)
					result[result_idx++] = false;
			}
		}

		return result;
	}*/

	string to_string() const
	{
		string res;

		for (size_t i = 0; i < p.size(); ++i)
		{
			if (p[i])
			{
				if (i)
					res = "x^" + std::to_string(i) + " + " + res;
				else
					res = "1";
			}
		}

		return res;
	}

	uint32_t get_msbit() const
	{
		uint32_t idx = 0;
		for (size_t i = 0; i < p.size(); ++i)
		{
			if (p[i])
				idx = i;
		}
		return idx+1;
	}
};

int main()
{
	string filename;
	int t;

	cout << "Filename: ";
	cin >> filename;

	cout << "Expected error rate: ";
	cin >> t;

	if (t < 1 || t > 8)
	{
		cout << "EER can be from 1 to 8" << endl;
		return 1;
	}

	int c = 2 * t;
	int e = pow(2, t) - 1;
	while (isPrime(c) || c % e == 0)
		++c;

	int n = NOK(e, c); // длина кода Файра
	int r = t + c; // длина избыточной информации
	int payload_len = n - r; // длина полезной информации

	// старшая степень не меньше длины исправляемого пакета (t)
	MyPoly gx;
	gx.p = gxes[t][0];

	MyPoly xc;
	xc.p[0] = 1;
	xc.p[c] = 1;

	MyPoly gxf = gx * xc;

	MyPoly msb_poly;
	msb_poly.p[gxf.get_msbit() - 1] = 1;
	
	MyPoly info;

	ifstream file(filename, ios::binary);
	ofstream output_file(filename + ".out", ios::binary);
	char in_byte = 0, out_byte = 0;
	int bits_read = 0, bits_wrote = 0;
	while (file.get(in_byte))
	{
		for (int i = 7; i >= 0; i--)
		{
			bool bit = (in_byte >> i) & 1;
			info.p |= bit;
			bits_read++;
			if (bits_read == payload_len)
			{
				MyPoly qx = info * msb_poly;
				MyPoly cx = (qx % gxf) + qx;

				for (int j = n - 1; j >= 0; j--)
				{
					// write bit
					out_byte |= cx.p[j];
					bits_wrote++;

					if (bits_wrote == 8)
					{
						output_file.put(out_byte);
						out_byte = 0;
						bits_wrote = 0;
					}

					out_byte <<= 1;
				}

				bits_read = 0;
				info.p = 0;
			}
			info.p <<= 1;
		}
	}

	/*size_t buf_size = 8 * (n - r);
	vector<uint8_t> buffer(buf_size);
	for (streamsize offset = 0; offset < filesize; offset += buf_size)
	{
		size_t this_buf_size = min((size_t) filesize - offset, buf_size);
		file.read((char*) buffer.data(), this_buf_size);

		size_t bit_offset = 0;
		for (uint8_t byte : buffer)
		{

		}
	}*/
}
