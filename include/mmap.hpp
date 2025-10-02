#ifndef CPP_MMAP_HPP_1753109898
#define CPP_MMAP_HPP_1753109898

#include <system_error>
#include <chrono>
#include <filesystem>
#include <string>
#include <exception>
#include <sys/file.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdint>
#include <unistd.h>
#include <sys/mman.h>

namespace fs = std::filesystem;

namespace mymm {

    inline void set_error_code_or_throw(std::error_code* ec, const std::string& msg = "")
    {
      if (ec)
      {
        ec->assign(errno, std::system_category());
      }
      else
      {
        std::error_code code;
        code.assign(errno, std::system_category());
        throw std::system_error(code, msg);
      }
    }

    using fd_handle_type = int;
    constexpr fd_handle_type invalid_fd = -1;

    enum class open_mode
    {
      read,   // r
      write,  // w
      append, // a
      read_binary, // rb
      write_binary, // wb
      append_binary, // ab
      read_write, // r+
      write_read, // w+
      append_rw, // a+
      read_write_binary, // rb+
      write_read_binary, // wb+
      append_rw_binary // ab+
    };

    inline fd_handle_type open_fd(fs::path path, int open_flags, int mode = 0644, std::error_code* ec = nullptr)
    {
      fd_handle_type handle = ::open(path.c_str(), open_flags, mode);
      if (handle == invalid_fd)
      {
        set_error_code_or_throw(ec, "open_fd() -> '" + path.string() + "'");
      }
      return handle;
    }
    
    inline fd_handle_type open_fd(fs::path path, open_mode om, int mode = 0644, std::error_code* ec = nullptr)
    {
      int flags = (om == open_mode::read) ? O_RDONLY : (O_RDWR | O_CREAT);
      return open_fd(path, flags, mode, ec);
    }

    inline void close_fd(fd_handle_type fd, std::error_code* ec = nullptr)
    {
      int res = ::close(fd);
      if (res == -1)
      {
        set_error_code_or_throw(ec, "close_fd() -> " + std::to_string(fd));
      }
    }

    inline std::size_t fd_count()
    {
      return std::distance(std::filesystem::directory_iterator("/proc/self/fd"), std::filesystem::directory_iterator{});
    }

    class fdw
{
public:
    fdw() : m_handle(invalid_fd), m_own(false) {}
    fdw(const fdw&) = delete;
    fdw& operator=(const fdw&) = delete;

    // Move constructor — transfers ownership
    fdw(fdw&& other) noexcept
        : m_handle(other.m_handle), m_own(other.m_own)
    {
        other.m_handle = invalid_fd;
        other.m_own = false;
    }

    // Move assignment — transfers ownership safely
    fdw& operator=(fdw&& other) noexcept
    {
        if (this != &other)
        {
            close(); // close current fd if owned
            m_handle = other.m_handle;
            m_own = other.m_own;
            other.m_handle = invalid_fd;
            other.m_own = false;
        }
        return *this;
    }

    fdw(fd_handle_type fd, bool own = true)
        : m_handle(fd), m_own(own)
    {}

    void set(fd_handle_type fd, bool own = true)
    {
        m_handle = fd;
        m_own = own;
    }

    ~fdw()
    {
        if (is_open())
            close();
    }

    operator bool() const noexcept
    {
        return is_open();
    }

    [[nodiscard]] fd_handle_type value() const noexcept
    {
        return m_handle;
    }

    void release() noexcept
    {
        m_own = false;
    }

    void acquire() noexcept
    {
        m_own = true;
    }

    [[nodiscard]] bool is_open() const noexcept
    {
        return m_handle != invalid_fd;
    }

    void close(bool force = false)
    {
        if ((m_own || force) && is_open())
        {
            close_fd(m_handle);
            m_handle = invalid_fd;
        }
    }

    private:
        fd_handle_type m_handle {invalid_fd};
        bool m_own {true};
    };

    [[nodiscard]] inline std::size_t get_page_size()
    {
      static std::size_t page_size = sysconf(_SC_PAGE_SIZE);
      return page_size;
    }

    [[nodiscard]] inline std::size_t make_aligned(std::size_t offset) noexcept
    {
      return offset & ~(get_page_size() - 1);
    }

    [[nodiscard]] inline std::size_t make_aligned_up(std::size_t offset) noexcept
    {
      return make_aligned(offset + get_page_size() - 1);

    }

   class filestat
    {
      public:
        filestat() = default;

        filestat(const std::string& path)
        {
          if (::stat(path.c_str(), &m_stats) == -1)
            set_error_code_or_throw(nullptr, "filestat::open() -> '" + path + "'");
        }

        filestat(int fd)
        {
          if (::fstat(fd, &m_stats) == -1)
            set_error_code_or_throw(nullptr, "filestat::open() -> fd " + std::to_string(fd));
        }

        bool open(const std::string& path)
        {
          return ::stat(path.c_str(), &m_stats) != -1;
        }

        bool open(int fd)
        {
          return ::fstat(fd, &m_stats) != -1;
        }

        [[nodiscard]] inline dev_t dev() const noexcept
        {
          return m_stats.st_dev;
        }

        [[nodiscard]] inline ino_t ino() const noexcept
        {
          return m_stats.st_ino;
        }

        [[nodiscard]] inline mode_t mode() const noexcept
        {
          return m_stats.st_mode;
        }

        [[nodiscard]] inline nlink_t nlink() const noexcept
        {
          return m_stats.st_nlink;
        }

        [[nodiscard]] inline uid_t uid() const noexcept
        {
          return m_stats.st_uid;
        }

        [[nodiscard]] inline gid_t gid() const noexcept
        {
          return m_stats.st_gid;
        }

        [[nodiscard]] inline dev_t rdev() const noexcept
        {
          return m_stats.st_rdev;
        }

        [[nodiscard]] inline off_t size() const noexcept
        {
          return m_stats.st_size;
        }

        [[nodiscard]] inline blksize_t blksize() const noexcept
        {
          return m_stats.st_blksize;
        }

        [[nodiscard]] inline blkcnt_t blocks() const noexcept
        {
          return m_stats.st_blocks;
        }

        template<typename Duration = std::chrono::seconds>
        [[nodiscard]] inline Duration access_time() const noexcept
        {
          return std::chrono::duration_cast<Duration>(
            std::chrono::nanoseconds(m_stats.st_atim.tv_nsec)
          );
        }

        template<typename Duration = std::chrono::seconds>
        [[nodiscard]] inline Duration modif_time() const noexcept
        {
          return std::chrono::duration_cast<Duration>(
            std::chrono::nanoseconds(m_stats.st_mtim.tv_nsec)
          );
        }

        template<typename Duration = std::chrono::seconds>
        [[nodiscard]] inline Duration cs_time() const noexcept
        {
          return std::chrono::duration_cast<Duration>(
            std::chrono::nanoseconds(m_stats.st_ctim.tv_nsec)
          );
        }

      private:
        struct stat m_stats;
    };

    struct mem_mapping
    {
      char* data {nullptr};
      std::size_t size {0};
      std::size_t mapped_size {0};
    };

    enum class mmap_flag : int
{
    shared = MAP_SHARED,
    priv = MAP_PRIVATE,
    
#ifdef MAP_SHARED_VALIDATE
    shared_validate = MAP_SHARED_VALIDATE,
#else
    shared_validate = MAP_SHARED, // fallback
#endif

#ifdef MAP_ANON
    anon = MAP_ANON,
#else
    anon = MAP_ANONYMOUS, // fallback if MAP_ANON not defined
#endif

    anonymous = MAP_ANONYMOUS,

#ifdef MAP_DENYWRITE
    denywrite = MAP_DENYWRITE,
#endif

#ifdef MAP_EXECUTABLE
    executable = MAP_EXECUTABLE,
#endif

    file = MAP_FILE,

#ifdef MAP_FIXED
    fixed = MAP_FIXED,
#else
    fixed = MAP_SHARED, // fallback
#endif

#ifdef MAP_FIXED_NOREPLACE
    fixed_noreplace = MAP_FIXED_NOREPLACE,
#else
    fixed_noreplace = MAP_FIXED, // fallback
#endif

#ifdef MAP_GROWSDOWN
    growsdown = MAP_GROWSDOWN,
#endif

#ifdef MAP_HUGETLB
    hugetlb = MAP_HUGETLB,
#endif

#ifdef MAP_LOCKED
    locked = MAP_LOCKED,
#endif

#ifdef MAP_NONBLOCK
    nonblock = MAP_NONBLOCK,
#endif

#ifdef MAP_NORESERVE
    noreserve = MAP_NORESERVE,
#endif

#ifdef MAP_POPULATE
    populate = MAP_POPULATE,
#endif

#ifdef MAP_STACK
    stack = MAP_STACK,
#endif

#ifdef MAP_SYNC
    sync = MAP_SYNC,
#endif

#ifdef MAP_HUGE_SHIFT
    huge_shift = MAP_HUGE_SHIFT,
#endif
};

    enum class mmap_prot : int
    {
      exec = PROT_EXEC,
      read = PROT_READ,
      write = PROT_WRITE,
      none = PROT_NONE,
      read_write = PROT_READ | PROT_WRITE,
    };

    enum class mmap_sync : int
    {
      async = MS_ASYNC,
      sync = MS_SYNC,
      invalidate = MS_INVALIDATE,
      sync_invalidate = MS_SYNC | MS_INVALIDATE,
      async_invalidate = MS_ASYNC | MS_INVALIDATE
    };

    enum class mmap_advise : int
    {
      normal = POSIX_MADV_NORMAL,
      sequential = POSIX_MADV_SEQUENTIAL,
      random = POSIX_MADV_RANDOM,
      willneed = POSIX_MADV_WILLNEED,
      dontneed = POSIX_MADV_DONTNEED,
    };

    inline mem_mapping make_mem_mapping(fd_handle_type handle, mmap_prot mp, mmap_flag mf, std::size_t size, std::size_t offset, std::error_code* ec = nullptr)
    {
      if (size == 0)
      {
        size = filestat(handle).size();
      }

      std::size_t align = make_aligned(offset);
      std::size_t map_size = (align + size) - offset;

      char* ptr = static_cast<char*>(::mmap(
        nullptr, map_size, static_cast<int>(mp), static_cast<int>(mf), handle, align
      ));

      if (ptr == MAP_FAILED)
      {
        set_error_code_or_throw(ec, "make_mem_mapping(): unable to map fd '" + std::to_string(handle) + "'");
        return mem_mapping { nullptr, 0, 0 };
      }

      return mem_mapping { ptr + offset - align, size, map_size };
    }

    inline void close_mem_mapping(mem_mapping* mapping, std::error_code* ec = nullptr)
    {
      char* start = mapping->data - (mapping->size - mapping->mapped_size);

      int res = ::munmap(start, mapping->size);

      if (res == -1)
      {
        set_error_code_or_throw(ec, "close_mem_mapping(): 'munmap' failed.");
      }

      mapping->data = nullptr;
    }

    inline void sync_mem_mapping(mem_mapping* mapping, mmap_sync ms, std::error_code* ec = nullptr)
    {
      char* start = mapping->data - (mapping->mapped_size - mapping->size);
      int res = ::msync(start, mapping->mapped_size, static_cast<int>(ms));

      if (res == -1)
      {
        set_error_code_or_throw(ec, "sync_mem_mapping(): 'msync' failed.");
      }
    }

    template<typename T, enum open_mode opm>
    class mmap
    {
      public:
        using value_type = T;
        using size_type = std::size_t;
        using reference = value_type&;
        using const_reference = const value_type&;
        using pointer = value_type*;
        using const_pointer = const value_type*;

        using iterator = pointer;
        using const_iterator = const_pointer;

      public:
        mmap(const mmap&) = delete;
        mmap& operator=(const mmap&) = delete;
        mmap(mmap&&) = default;
        mmap& operator=(mmap&&) = default;

        mmap() = default;

        mmap(const fs::path& path, std::size_t size, std::size_t offset, mmap_flag flag = mmap_flag::shared)
        {
          open(path, size, offset, flag);
        }

        mmap(fd_handle_type fd, std::size_t size, std::size_t offset, mmap_flag flag = mmap_flag::shared, bool own = false)
        {
          open(fd, size, offset, flag);
        }

        void open(const fs::path& path, std::size_t size, std::size_t offset, mmap_flag flag = mmap_flag::shared)
        {
          m_fdw.set(open_fd(path, opm));
          init(flag, size, offset);
        }

        void open(fd_handle_type fd, std::size_t size, std::size_t offset, mmap_flag flag = mmap_flag::shared, bool own = false)
        {
          m_fdw.set(fd, own);
          init(flag, size, offset);
        }

        template<open_mode om_ = opm, typename = std::enable_if_t<om_ == open_mode::write, void>>
        void create(const fs::path& path, std::size_t size, std::size_t offset, mmap_flag flag = mmap_flag::shared)
        {
          m_fdw.set(open_fd(path, opm));
          init(flag, size, offset, true);
        }

        void close_mapping()
        {
          if (is_mapped())
          {
            close_mem_mapping(&m_mapping);
          }
        }

        void close_fd()
        {
          if (is_open())
          {
            m_fdw.close();
          }
        }

        void close()
        {
          close_mapping();
          close_fd();
        }

        void advise(mmap_advise ma, std::size_t offset = 0, std::size_t length = 0, std::error_code* ec = nullptr)
        {
          if (length == 0)
          {
            length = mapped_size();
          }

          if (offset + length > mapped_size())
          {
            length = mapped_size() - offset;
          }

          int res = ::posix_madvise(m_mapping.data + offset * sizeof(value_type), length * sizeof(value_type), static_cast<int>(ma));

          if (res != 0)
          {
            set_error_code_or_throw(ec, "posix_madvise() failed");
          }
        }

      private:
        void init(mmap_flag flag, std::size_t size, std::size_t offset, bool trunc = false)
        {
          size = size * sizeof(value_type);
          offset = offset * sizeof(value_type);

          if (trunc)
        {
            int res = ::ftruncate(m_fdw.value(), size);
            if (res != 0) {
                set_error_code_or_throw(nullptr, "mmap::init()");
            }
        }

          if constexpr(opm == open_mode::read)
          {
            m_mapping = make_mem_mapping(m_fdw.value(), mmap_prot::read, flag, size, offset);
          }
          else
          {
            m_mapping = make_mem_mapping(m_fdw.value(), mmap_prot::read_write, flag, size, offset);
          }
        }

      public:

        void sync(mmap_sync ms = mmap_sync::sync)
        {
          sync_mem_mapping(&m_mapping, ms);
        }

        [[nodiscard]] size_type size() const noexcept
        {
          return m_mapping.size / sizeof(value_type);
        }

        [[nodiscard]] size_type size_byte() const noexcept
        {
          return m_mapping.size;
        }

        [[nodiscard]] size_type mapped_size() const noexcept
        {
          return m_mapping.mapped_size / sizeof(value_type);
        }

        [[nodiscard]] size_type mapped_size_byte() const noexcept
        {
          return m_mapping.mapped_size;
        }

        [[nodiscard]] bool is_mapped() const noexcept
        {
          return m_mapping.data != nullptr;
        }

        [[nodiscard]] bool is_open() const noexcept
        {
          return m_fdw.is_open();
        }

        [[nodiscard]] const char* bytes() const noexcept
        {
          return m_mapping.data;
        }

        [[nodiscard]] const_pointer data() const noexcept
        {
          return reinterpret_cast<const_pointer>(bytes());
        }

        [[nodiscard]] const_iterator begin() const noexcept
        {
          return data();
        }

        [[nodiscard]] const_iterator end() const noexcept
        {
          return data() + mapped_size();
        }

        [[nodiscard]] const_iterator cbegin() const noexcept
        {
          return data();
        }

        [[nodiscard]] const_iterator cend() const noexcept
        {
          return data() + mapped_size();
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] char* bytes() noexcept
        {
          return m_mapping.data;
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] pointer data() noexcept
        {
          return reinterpret_cast<pointer>(bytes());
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] iterator begin() noexcept
        {
          return data();
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] iterator end() noexcept
        {
          return data() + mapped_size();
        }

        [[nodiscard]] const_reference operator[](const size_type index) const noexcept
        {
          return data()[index];
        }

        [[nodiscard]] const_reference at(const size_type index) const
        {
          if (index >= mapped_size())
            throw std::out_of_range("mmap::at()");
          return (*this)[index];
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] reference operator[](const size_type index) noexcept
        {
          return data()[index];
        }

        template<open_mode opm_ = opm, typename = std::enable_if_t<opm_ == open_mode::write>>
        [[nodiscard]] reference at(const size_type index)
        {
          if (index >= mapped_size())
            throw std::out_of_range("mmap::at()");
          return (*this)[index];
        }

        [[nodiscard]] const unsigned char& at_byte(const size_type index) const noexcept
        {
          return static_cast<const unsigned char*>(data())[index];
        }

        [[nodiscard]] unsigned char& at_byte(const size_type index)
        {
          if (index >= mapped_size_byte())
            throw std::out_of_range("mmap::at_byte()");
          return static_cast<unsigned char*>(data())[index];
        }

        template<template<typename> typename Container>
        void collect(Container<T>& ct) const
        {
          ct.insert(ct.end(), begin(), end());
        }

      private:
        fdw m_fdw {-1, false};
        mem_mapping m_mapping;
    };

    template<typename T>
    using immap = mmap<T, open_mode::read>;

    template<typename T>
    using ommap = mmap<T, open_mode::write>;

}


#endif /* end of include guard: CPP_MMAP_HPP_1753109898 */
