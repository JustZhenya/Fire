#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <mpi.h>
#include <thread>

using namespace std;

// Наибольший общий делитель
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

class VecPoly
{
private:
	static vector<bool> my_xor(const vector<bool>& a, const vector<bool>& b, size_t offset = 1)
	{
		vector<bool> res;

		for (size_t i = offset; i < min(a.size(), b.size()); i++)
			res.push_back(a[i] != b[i]);

		return res;
	}

public:

	vector<bool> vec;

	VecPoly(size_t len) : vec(len, false) { }
	VecPoly(const vector<bool>& _vec) : vec(_vec) { }

	operator bool() const
	{
		for (bool b : vec)
			if (b)
				return true;

		return false;
	}

	VecPoly operator>>(int offset) const
	{
		vector<bool> res(vec.size(), 0);
		if (vec.size() >= offset)
		{
			copy(vec.begin() + offset, vec.end(), res.begin());
		}
		return res;
	}

	VecPoly operator<<(int offset) const
	{
		vector<bool> res(vec.size(), 0);
		if (vec.size() >= offset)
		{
			copy(vec.begin(), vec.end() - offset, res.begin() + offset);
		}
		return res;
	}

	VecPoly operator&(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < min(vec.size(), other.vec.size()); i++)
			res.push_back(vec[i] && other.vec[i]);

		return res;
	}

	VecPoly operator|(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < min(vec.size(), other.vec.size()); i++)
			res.push_back(vec[i] || other.vec[i]);

		return res;
	}

	VecPoly operator^(const VecPoly& other) const
	{
		return my_xor(vec, other.vec, 0);
	}

	VecPoly operator*(const VecPoly& other) const
	{
		vector<bool> res;

		for (size_t i = 0; i < other.vec.size(); ++i)
		{
			for (size_t j = 0; j < vec.size(); ++j)
			{
				if (other.vec[i] && vec[j])
				{
					size_t offset = i + j;
					if (res.size() < offset + 1)
						res.resize(offset + 1, false);
					res[offset] = true;
				}
			}
		}

		return res;
	}

	VecPoly operator%(const VecPoly& divisor) const
	{
		size_t pick = divisor.vec.size();

		VecPoly tmp(pick);
		copy(vec.begin(), vec.begin() + min(pick, vec.size()), tmp.vec.begin());

		//cout << "begin: " << tmp.to_string() << endl;

		for (; pick < vec.size(); pick++)
		{
			//cout << tmp.to_string() << " " << pick << endl;
			if (tmp.vec[0])
			{
				tmp.vec = my_xor(divisor.vec, tmp.vec);
			}
			else
			{
				VecPoly t(pick);
				tmp.vec = my_xor(t.vec, tmp.vec);
			}
			tmp.vec.push_back(vec[pick]);
			//cout << tmp.to_string() << " " << pick+1 << endl;
		}

		if (tmp.vec[0])
		{
			tmp.vec = my_xor(divisor.vec, tmp.vec);
		}
		else
		{
			VecPoly t(pick);
			tmp.vec = my_xor(t.vec, tmp.vec);
		}

		return tmp;
	}

	string to_string() const
	{
		string res;

		for (bool b : vec)
			res += b ? "1" : "0";

		return res;
	}

	void crshift(int len)
	{
		*this = (*this >> len) | (*this << (int)(vec.size() - len));
	}

	void clshift(int len)
	{
		*this = (*this << len) | (*this >> (int)(vec.size() - len));
	}
};

class FireEncDec
{
public:
	// неприводимые полиномы g(x)
	static const map<int, vector<vector<bool>>> gxes;

	int n; // блоковая длина кода Файра
	int r; // длина проверочной группы в кодовой комбинации
	int payload_len; // длина полезной информации
	size_t gxf_msb; // старшая степень полинома g(x)f

	// полиномы для декодирования
	VecPoly gx;
	VecPoly xc;

	// полином для кодирования
	VecPoly gxf;

public:
	// P - длина исправляемого пакета
	FireEncDec(int P) : n(), r(), payload_len(), gxf_msb(), gx(0), xc(0), gxf(0)
	{
		int t = P; // старшая степень t не меньше длины исправляемого пакета P
		int c = 2 * P; // Старшая степень с второго полинома должна равняться или превышать удвоенную длину предполагаемого пакета ошибок Р
		int e = pow(2, t) - 1; // При этом с не должно делиться нацело на число е, где е = 2^t - 1
		while (c % e == 0)
			c++;

		gx = gxes.at(t)[0];

		xc.vec.resize(c + 1, false);
		xc.vec[0] = 1;
		xc.vec[c] = 1;

		gxf = gx * xc; // образующий полином кода Файра G(X)Ф определяется произведением двух примитивных полиномов: G(X)Ф = G(X)(Xс+1)
		gxf_msb = gxf.vec.size() - 1; // старшая степень полинома g(x)f

		n = e * c; // блоковая длина кода Файра будет являться произведением чисел e и с: n = е*с
		r = t + c; // длина проверочной группы r в кодовой комбинации будет являться суммой старших степеней обоих образующих полиномов: r = t+c
		payload_len = n - r; // длина полезной информации

		/*cout << "P = " << P << endl;
		cout << "t = " << t << endl;
		cout << "c = " << c << endl;
		cout << "e = " << e << endl;
		cout << "gx = " << gx.to_string() << endl;
		cout << "xc = " << xc.to_string() << endl;
		cout << "gxf = " << gxf.to_string() << endl;
		cout << "gxf_msb = " << gxf_msb << endl;
		cout << "n = " << n << endl;
		cout << "r = " << r << endl;
		cout << "n-r = " << payload_len << endl;*/
	}

	static long filesize(const string& filename)
	{
		struct stat stat_buf;
		int rc = stat(filename.c_str(), &stat_buf);
		return rc == 0 ? stat_buf.st_size : -1;
	}

	void encode_block(VecPoly& info) const
	{
		info.vec.resize(info.vec.size() + gxf_msb, 0); // Производим увеличение степени информационного полинома I(X) на старшую степень образующего полинома G(X)Ф
		VecPoly remainder = info % gxf; // Выполняем деление расширенного информационного полинома Q(X) на образующий полином G(X)Ф
		info.vec.erase(info.vec.begin() + payload_len, info.vec.end());
		info.vec.insert(info.vec.end(), remainder.vec.begin(), remainder.vec.end());
	}

	void decode_block(VecPoly& cx) const
	{
		// Декодирование осуществляется путем деления кодового полинома С(Х) последовательно на обе компоненты образующего полинома кода Файра G(X)Ф: на (Х2+Х+1) и на (Х6+1).
		VecPoly info_remainder1 = cx % gx;
		VecPoly info_remainder2 = cx % xc;

		if (info_remainder1 && info_remainder2)
		{
			// обнаружена ошибка
			cout << "Found error!" << endl;

			VecPoly cx_orig = cx;
			size_t shift_count = 0;
			while (info_remainder1.vec != info_remainder2.vec && shift_count < n)
			{
				cx.crshift(1);
				info_remainder1 = cx % gx;
				info_remainder2 = cx % xc;
				shift_count++;

				// удалить нули в начале вектора
				while (!info_remainder1.vec.empty() && !info_remainder1.vec[0])
					info_remainder1.vec.erase(info_remainder1.vec.begin());

				while (!info_remainder2.vec.empty() && !info_remainder2.vec[0])
					info_remainder2.vec.erase(info_remainder2.vec.begin());
			}

			info_remainder1.vec.insert(info_remainder1.vec.begin(), cx_orig.vec.size() - info_remainder1.vec.size(), 0);
			info_remainder1.clshift(shift_count);
			cx = cx_orig ^ info_remainder1;
		}

		cx.vec.erase(cx.vec.end() - r, cx.vec.end()); // отбросить проверочный код и поделить на x^6 (на самом деле x^4)
	}
};

const map<int, vector<vector<bool>>> FireEncDec::gxes = {
	{2, {{1, 1, 1}}},
	{3, {{1, 0, 1, 1}}},
	{4, {{1, 0, 0, 1, 1}, {1, 1, 1, 1, 1}}},
	{5, {{1, 0, 0, 1, 0, 1}}},
	{6, {{1, 0, 0, 0, 0, 1, 1}, {1, 0, 0, 1, 0, 0, 1}, {1, 1, 0, 1, 1, 1, 1}}},
	{7, {{1, 1, 1, 0, 0, 1, 1, 1}, {1, 0, 0, 0, 0, 0, 1, 1}, {1, 0, 0, 1, 1, 1, 0, 1}}}
};

class MPIController
{
private:
	const FireEncDec& fire;
	string filename_in;
	string filename_out;
	size_t payload_total_bytes;
	size_t out_bytes;
	int mode;

public:
	MPIController(const FireEncDec& _fire, const string& _filename_in, const string& _filename_out, size_t _payload_total_bytes, int _mode) :
		fire(_fire),
		filename_in(_filename_in),
		filename_out(_filename_out),
		payload_total_bytes(_payload_total_bytes),
		out_bytes(0),
		mode(_mode)
	{
		 
	}

	void send_data()
	{
		ifstream in_file(filename_in, ios::binary);

		if (mode == 1)
			in_file.read((char*)&out_bytes, sizeof(out_bytes));

		size_t in_block_size = mode == 0 ? fire.payload_len : fire.n; // in bits
		size_t out_block_size = mode == 0 ? fire.n : fire.payload_len; // in bits
		size_t in_chunk_size = NOK(16, NOK(in_block_size, out_block_size)); // in bits
		size_t blocks_per_chunk = in_chunk_size / in_block_size;
		size_t out_chunk_size = out_block_size * blocks_per_chunk; // in bits
		size_t total_chunks = ceil(payload_total_bytes * 8.0 / in_chunk_size);

		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);

		size_t chunks_per_thread = ceil((double) total_chunks / (world_size - 1));

		cout << in_block_size << " " << out_block_size << " " << in_chunk_size << " " << chunks_per_thread << endl;

		size_t chunk_offset = 0;
		for (int thread_id = 1; thread_id < world_size; thread_id++)
		{
			size_t this_thread_chunks = min(chunks_per_thread, total_chunks - chunk_offset);

			vector<uint8_t> buf((in_chunk_size / 8) * this_thread_chunks);
			in_file.read((char*) buf.data(), buf.size());

			MPI_Send(&this_thread_chunks, 1, MPI_LONG_LONG, thread_id, 0, MPI_COMM_WORLD);
			MPI_Send(&chunk_offset, 1, MPI_LONG_LONG, thread_id, 0, MPI_COMM_WORLD);
			MPI_Send(buf.data(), buf.size(), MPI_UNSIGNED_CHAR, thread_id, 0, MPI_COMM_WORLD);

			chunk_offset += this_thread_chunks;
		}
	}

	void recv_data()
	{
		ofstream out_file(filename_out, ios::binary);

		if (mode == 0)
			out_file.write((char*)&payload_total_bytes, sizeof(payload_total_bytes));

		size_t in_block_size = mode == 0 ? fire.payload_len : fire.n; // in bits
		size_t out_block_size = mode == 0 ? fire.n : fire.payload_len; // in bits
		size_t in_chunk_size = NOK(16, NOK(in_block_size, out_block_size)); // in bits
		size_t blocks_per_chunk = in_chunk_size / in_block_size;
		size_t out_chunk_size = out_block_size * blocks_per_chunk; // in bits

		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);

		size_t bytes_written = 0;
		for (int thread_id = 1; thread_id < world_size; thread_id++)
		{
			size_t this_thread_chunks, chunk_offset;
			MPI_Recv(&this_thread_chunks, 1, MPI_LONG_LONG, thread_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&chunk_offset, 1, MPI_LONG_LONG, thread_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			vector<char> buf((out_chunk_size / 8) * this_thread_chunks);
			MPI_Recv(buf.data(), buf.size(), MPI_UNSIGNED_CHAR, thread_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			size_t bytes_to_write = mode == 0 ? buf.size() : min(buf.size(), out_bytes - bytes_written);
			out_file.write(buf.data(), bytes_to_write);
			bytes_written += bytes_to_write;
		}
	}

	void run()
	{
		thread send_thr(&MPIController::send_data, this); // send chunks to all threads
		thread recv_thr(&MPIController::recv_data, this); // recv chunks from threads
		send_thr.join();
		recv_thr.join();
	}
};

class MPIWorker
{
private:
	const FireEncDec& fire;
	int mode;

public:
	MPIWorker(const FireEncDec& _fire, int _mode) :
		fire(_fire),
		mode(_mode)
	{

	}

	void run()
	{
		size_t in_block_size = mode == 0 ? fire.payload_len : fire.n; // in bits
		size_t out_block_size = mode == 0 ? fire.n : fire.payload_len; // in bits
		size_t in_chunk_size = NOK(16, NOK(in_block_size, out_block_size)); // in bits

		size_t this_thread_chunks, chunk_offset;
		MPI_Recv(&this_thread_chunks, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&chunk_offset, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		vector<uint8_t> buf((in_chunk_size / 8) * this_thread_chunks);
		MPI_Recv(buf.data(), buf.size(), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		vector<bool> bitchunk;
		bitchunk.reserve(in_chunk_size);
		for (uint8_t byte_data : buf)
		{
			for (int j = 7; j >= 0; j--)
				bitchunk.push_back((byte_data >> j) & 1);
		}

		buf.clear();
		//buf.reserve();

		int bits_wrote = 0;
		uint8_t out_byte = 0;
		size_t blocks = bitchunk.size() / in_block_size;
		for (size_t i = 0; i < blocks; i++)
		{
			VecPoly info(in_block_size);
			size_t in_bit_offset = i * in_block_size;
			copy(bitchunk.begin() + in_bit_offset, bitchunk.begin() + in_bit_offset + in_block_size, info.vec.begin());

			if (mode == 0)
				fire.encode_block(info);
			else
				fire.decode_block(info);

			for (bool b : info.vec)
			{
				out_byte |= b ? 1 : 0;
				bits_wrote++;

				if (bits_wrote == 8)
				{
					buf.push_back(out_byte);
					out_byte = 0;
					bits_wrote = 0;
				}

				out_byte <<= 1;
			}
		}

		MPI_Send(&this_thread_chunks, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&chunk_offset, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
		MPI_Send(buf.data(), buf.size(), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	}
};

int main(int argc, char** argv)
{
	if (MPI_Init(&argc, &argv))
	{
		cerr << "Failed to init MPI" << endl;
		return 1;
	}

	string filename = argv[1];
	int P = stoi(argv[2]);
	int mode = stoi(argv[3]);

	if (mode == 0)
	{
		cout << "Encoding..." << flush;
	}
	else if (mode == 1)
	{
		cout << "Decoding..." << flush;
	}
	else
	{
		cerr << "Invalid Mode" << endl;
		MPI_Finalize();
		return 1;
	}

	FireEncDec fire(P);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		size_t bytes = FireEncDec::filesize(filename);
		if(mode == 1)
			bytes -= sizeof(size_t);
		string out_filename = filename + (mode == 0 ? ".enc" : ".dec");
		MPIController mpic(fire, filename, out_filename, bytes, mode);
		mpic.run();
	}
	else
	{
		MPIWorker mpiw(fire, mode);
		mpiw.run();
	}

	cout << rank << " done" << endl;

	MPI_Finalize();
}
