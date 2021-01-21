
#pragma once


namespace MyOptional
{
	const int CODE_CORRECT = 0;
	const int CODE_BAD_SIZE = 1;
	const int CODE_BAD_VALUE = 2;
	const int CODE_SINGULAR_MATRIX = 3;
	const int CODE_FILE_NOT_FOUND = 4;
	
	template<typename T>
	class COptional
	{
	public:
		COptional<T>(T* item, int code = CODE_CORRECT);
		COptional<T>(int code);
		COptional<T>();

		T* get_value();
		int get_code();

		bool is_correct();

		void set_value(T* value);
		void set_code(int code);
	private:
		T* _value;
		int _code;
	};

	template <typename T>
	COptional<T>::COptional(T* item, const int code): _value(item), _code(code)
	{
	}

	template <typename T>
	COptional<T>::COptional(const int code): _value(nullptr), _code(code)
	{
	}

	template <typename T>
	COptional<T>::COptional():_value(nullptr), _code(CODE_CORRECT)
	{
	}

	template <typename T>
	T* COptional<T>::get_value()
	{
		return _value;
	}

	template <typename T>
	int COptional<T>::get_code()
	{
		return _code;
	}

	template <typename T>
	bool COptional<T>::is_correct()
	{
		return _code == CODE_CORRECT;
	}

	template <typename T>
	void COptional<T>::set_value(T* value)
	{
		_value = value;
	}

	template <typename T>
	void COptional<T>::set_code(int code)
	{
		_code = code;
	}
}
