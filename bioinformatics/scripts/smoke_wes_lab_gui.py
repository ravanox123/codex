"""Exercise real Tk screens against the saved sample; optionally capture the app window."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from wes_lab.app import App
from wes_lab.data import Project

def main():
    root=Path(__file__).resolve().parents[1]
    app=App(restore=False);app.loaded(Project(root/'results/wes_lab/571566'));errors=[]
    app.report_callback_exception=lambda *args:errors.append(args)
    try:
        app.update()
        for page in ('Projects','Files','Overview','Variants','Coverage','Reports','Run analysis'):
            app.show_page(page);app.update()
            if page=='Projects':assert app.people_table.get_children()
            if page=='Files':
                assert app.person.data['name']=='571566'
                assert len(app.files_tabs.tabs())==3
                assert any(r['name'].endswith('.bam') for r in app.file_entries.values())
                assert app.analysis_entries
            if page=='Variants':
                app.search_var.set('FDXR');app.group_var.set('Clinical assertion');app.apply_filters();app.update()
                assert len(app.table.get_children())==1
                iid=app.table.get_children()[0];app.table.selection_set(iid);app.variant_selected();assert 'p.Arg386Trp' in app.detail_text.get('1.0','end')
            if page=='Coverage':
                app.coverage_query.set('SYNGAP1');app.refresh_coverage();app.update();assert app.coverage_tree.get_children()==('SYNGAP1',)
            print('GUI screen OK:',page,flush=True)
        assert not errors,errors
        app.show_page('Files');app.update()
        if '--capture' in sys.argv:capture(app,root/'results/wes_lab/app-preview.png')
        print('GUI smoke passed',flush=True)
    finally:app.destroy()

def capture(app,path):
    import ctypes as c,ctypes.util,struct,zlib
    lib=c.CDLL(ctypes.util.find_library('X11'))
    class XImage(c.Structure):
        _fields_=[('width',c.c_int),('height',c.c_int),('xoffset',c.c_int),('format',c.c_int),('data',c.c_void_p),('byte_order',c.c_int),('bitmap_unit',c.c_int),('bitmap_bit_order',c.c_int),('bitmap_pad',c.c_int),('depth',c.c_int),('bytes_per_line',c.c_int),('bits_per_pixel',c.c_int),('red_mask',c.c_ulong),('green_mask',c.c_ulong),('blue_mask',c.c_ulong)]
    lib.XOpenDisplay.restype=c.c_void_p;lib.XOpenDisplay.argtypes=[c.c_char_p];display=lib.XOpenDisplay(None)
    if not display:raise RuntimeError('X11 display unavailable for preview capture.')
    lib.XGetImage.restype=c.POINTER(XImage);lib.XGetImage.argtypes=[c.c_void_p,c.c_ulong,c.c_int,c.c_int,c.c_uint,c.c_uint,c.c_ulong,c.c_int]
    lib.XGetPixel.argtypes=[c.POINTER(XImage),c.c_int,c.c_int];lib.XGetPixel.restype=c.c_ulong
    lib.XDestroyImage.argtypes=[c.POINTER(XImage)];lib.XCloseDisplay.argtypes=[c.c_void_p]
    w=app.winfo_width();h=app.winfo_height();im=lib.XGetImage(display,app.winfo_id(),0,0,w,h,c.c_ulong(-1).value,2)
    if not im:raise RuntimeError('Cannot capture the Tk window.')
    try:
        raw=bytearray()
        for y in range(h):
            raw.append(0)
            for x in range(w):
                p=lib.XGetPixel(im,x,y);raw.extend(((p>>16)&255,(p>>8)&255,p&255))
        def chunk(name,data):return struct.pack('!I',len(data))+name+data+struct.pack('!I',zlib.crc32(name+data)&0xffffffff)
        path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('!IIBBBBB',w,h,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(raw,6))+chunk(b'IEND',b''))
        print('Preview saved:',path,flush=True)
    finally:lib.XDestroyImage(im);lib.XCloseDisplay(display)

if __name__=='__main__':main()
