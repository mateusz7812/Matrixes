
#pragma once


namespace MyAlgebra
{
	const int CODE_CORRECT = 0;
	const int CODE_BAD_SIZE = 1;
	const int CODE_BAD_VALUE = 2;
	const int CODE_SINGULAR_MATRIX = 3;
	
	template<typename T>
	class Optional
	{
	public:
		Optional<T>(T* item, int code = CODE_CORRECT);
		Optional<T>(int code);
		Optional<T>();

		T* get_value();
		int get_code();

		void set_value(T* value);
		void set_code(int code);
	private:
		T* _value;
		int _code;
	};

	template <typename T>
	Optional<T>::Optional(T* item, const int code): _value(item), _code(code)
	{
	}

	template <typename T>
	Optional<T>::Optional(const int code): _value(nullptr), _code(code)
	{
	}

	template <typename T>
	Optional<T>::Optional():_value(nullptr), _code(CODE_CORRECT)
	{
	}

	template <typename T>
	T* Optional<T>::get_value()
	{
		return _value;
	}

	template <typename T>
	int Optional<T>::get_code()
	{
		return _code;
	}

	template <typename T>
	void Optional<T>::set_value(T* value)
	{
		_value = value;
	}

	template <typename T>
	void Optional<T>::set_code(int code)
	{
		_code = code;
	}
}
