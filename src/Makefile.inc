# 
# ------ Lists SPIDER source files for use by Makefile  ------
#        This is the definitive source listing for source files
#        compiled into SPIDER. Additional objects are also found
#        within the corresponding source files.
 
ELEMENTS =  \
 $(LIB)(spider.o)    \
 $(LIB)(openinln.o)  \
 $(LIB)(lunsethdr.o) \
 $(LIB)(lunsetmrchdr.o) \
 $(LIB)(fftw3.o)     \
 $(LIB)(lundoc.o)    \
 $(LIB)(apsh_tool.o) \
 $(LIB)(acrs_2.o)    \
 $(LIB)(acrs_2s.o)   \
 $(LIB)(acrs_3.o)    \
 $(LIB)(acrs_3s.o)   \
 $(LIB)(add.o)   \
 $(LIB)(ads.o)   \
 $(LIB)(addfac.o)\
 $(LIB)(adds.o)  \
 $(LIB)(adds_n.o)\
 $(LIB)(af.o)    \
 $(LIB)(afs.o)   \
 $(LIB)(akima.o) \
 $(LIB)(al10.o)  \
 $(LIB)(alphaint.o) \
 $(LIB)(alprbs.o)   \
 $(LIB)(alprbs_q.o) \
 $(LIB)(alrosf_q.o) \
 $(LIB)(alrosi_q.o) \
 $(LIB)(alrq.o)     \
 $(LIB)(alrqs.o)    \
 $(LIB)(alrq_m.o)   \
 $(LIB)(alrq_ms_quad.o) \
 $(LIB)(alrq_ms_fbs.o)  \
 $(LIB)(alrq_q.o) \
 $(LIB)(amoeba.o) \
 $(LIB)(ang.o)    \
 $(LIB)(ang_n.o)  \
 $(LIB)(aniso.o)  \
 $(LIB)(anisoe.o) \
 $(LIB)(anisof.o) \
 $(LIB)(anisof3.o)\
 $(LIB)(apcc.o)   \
 $(LIB)(apscc.o)  \
 $(LIB)(apshift.o)\
 $(LIB)(apshift_fou.o)\
 $(LIB)(ap_end.o)     \
 $(LIB)(ap_getangas.o)\
 $(LIB)(ap_getdata.o) \
 $(LIB)(apref_p.o)    \
 $(LIB)(apref_pm.o)   \
 $(LIB)(apsh_psc.o)   \
 $(LIB)(apsh_ps.o)    \
 $(LIB)(apsh_ss.o)    \
 $(LIB)(ap_stat.o)    \
 $(LIB)(apmaster.o)   \
 $(LIB)(apmaster_tool.o)\
 $(LIB)(aprings.o)      \
 $(LIB)(aprings_coef.o) \
 $(LIB)(aprings_trans.o)\
 $(LIB)(applyw.o)  \
 $(LIB)(applyws.o) \
 $(LIB)(arasq.o)   \
 $(LIB)(arith.o)   \
 $(LIB)(arithl.o)  \
 $(LIB)(arithsca.o)\
 $(LIB)(asta.o)    \
 $(LIB)(astcyl.o)  \
 $(LIB)(atmc.o)  \
 $(LIB)(atpk.o)  \
 $(LIB)(atsa.o)  \
 $(LIB)(atwn.o)  \
 $(LIB)(averg.o) \
 $(LIB)(averg3.o)\
 $(LIB)(bckc0.o) \
 $(LIB)(bckc2.o) \
 $(LIB)(bckpj.o) \
 $(LIB)(bcq.o)   \
 $(LIB)(besi0.o) \
 $(LIB)(besi1.o) \
 $(LIB)(betai.o) \
 $(LIB)(bfact.o) \
 $(LIB)(bigalloc.o) \
 $(LIB)(bigended.o) \
 $(LIB)(blank.o)  \
 $(LIB)(bldr.o)   \
 $(LIB)(blob.o)   \
 $(LIB)(blob_q.o) \
 $(LIB)(bmax_c.o) \
 $(LIB)(bmax_s.o) \
 $(LIB)(bmin_c.o) \
 $(LIB)(bmin_s.o) \
 $(LIB)(boxx.o)   \
 $(LIB)(bpcm.o)   \
 $(LIB)(bpcmp.o)  \
 $(LIB)(bpcq.o)   \
 $(LIB)(bpcqp.o)  \
 $(LIB)(bprp.o)   \
 $(LIB)(bprp3.o)  \
 $(LIB)(bprj2.o)  \
 $(LIB)(bps2.o)   \
 $(LIB)(bpwr.o)   \
 $(LIB)(bpwr_q.o) \
 $(LIB)(builds.o) \
 $(LIB)(buildm.o) \
 $(LIB)(buterp.o) \
 $(LIB)(calc.o)   \
 $(LIB)(cald.o)   \
 $(LIB)(cang.o)   \
 $(LIB)(ccc.o)    \
 $(LIB)(cconect.o)\
 $(LIB)(ccrd_2.o) \
 $(LIB)(ccrd_3.o) \
 $(LIB)(cnrd_2.o) \
 $(LIB)(cnrd_3.o) \
 $(LIB)(ccrs.o)   \
 $(LIB)(ccrs_2.o) \
 $(LIB)(ccrs_ph_3.o)\
 $(LIB)(cnrs_2.o)   \
 $(LIB)(cnrs_3.o)   \
 $(LIB)(cengr3.o)   \
 $(LIB)(cent_3ph.o) \
 $(LIB)(cent_ph.o)  \
 $(LIB)(cent_sym.o) \
 $(LIB)(cent_q.o)   \
 $(LIB)(cgr_3.o)    \
 $(LIB)(changeval.o) \
 $(LIB)(charinside.o)\
 $(LIB)(chava.o)   \
 $(LIB)(chisq.o)   \
 $(LIB)(chkinpq.o) \
 $(LIB)(chkstr.o)  \
 $(LIB)(chkstuff.o)\
 $(LIB)(clast.o)  \
 $(LIB)(clss.o)   \
 $(LIB)(cncalc.o) \
 $(LIB)(cnint3.o) \
 $(LIB)(cnscan.o) \
 $(LIB)(cnstuff.o)\
 $(LIB)(cntrce.o) \
 $(LIB)(cntrfl.o) \
 $(LIB)(cntur.o)  \
 $(LIB)(comp3dmad.o) \
 $(LIB)(conf.o)  \
 $(LIB)(conf1.o) \
 $(LIB)(conint.o)\
 $(LIB)(cop.o)   \
 $(LIB)(copy1.o) \
 $(LIB)(copybrix.o)    \
 $(LIB)(copyfrommrc.o) \
 $(LIB)(copytomrc.o)   \
 $(LIB)(copyemi.o) \
 $(LIB)(copyd.o)   \
 $(LIB)(copye.o)   \
 $(LIB)(copyf.o)   \
 $(LIB)(copyfromnik.o) \
 $(LIB)(copyfromtif.o) \
 $(LIB)(copyfromxplo.o)\
 $(LIB)(copypos.o)     \
 $(LIB)(copytojpg.o)   \
 $(LIB)(copytostk.o)   \
 $(LIB)(copytoxplor.o) \
 $(LIB)(copytovol.o)   \
 $(LIB)(copyr.o) \
 $(LIB)(copyu.o) \
 $(LIB)(corel.o) \
 $(LIB)(corr1.o) \
 $(LIB)(coupe.o) \
 $(LIB)(cppb.o)  \
 $(LIB)(cproj.o) \
 $(LIB)(cr180_2.o)   \
 $(LIB)(crcse1.o)    \
 $(LIB)(crcse2.o)    \
 $(LIB)(crcse3.o)    \
 $(LIB)(crosrmg.o)   \
 $(LIB)(crosrng.o)   \
 $(LIB)(crosrng_2.o) \
 $(LIB)(crosrng_trans.o) \
 $(LIB)(crosrng_ds.o)  \
 $(LIB)(crosrng_com.o) \
 $(LIB)(crosrng_ms.o)  \
 $(LIB)(crosrng_q.o)   \
 $(LIB)(crsm_2.o) \
 $(LIB)(cslice.o) \
 $(LIB)(csvms.o)  \
 $(LIB)(ctffind3.o)  \
 $(LIB)(cutxslice.o) \
 $(LIB)(cutyslice.o) \
 $(LIB)(cutzslice.o) \
 $(LIB)(date.o)    \
 $(LIB)(date_2k.o) \
 $(LIB)(d3dplt.o)  \
 $(LIB)(dchist.o)  \
 $(LIB)(decimate.o) \
 $(LIB)(dee.o)     \
 $(LIB)(defo001.o) \
 $(LIB)(defo003.o) \
 $(LIB)(defo1001.o) \
 $(LIB)(defo2001.o) \
 $(LIB)(defocus.o) \
 $(LIB)(deletf.o) \
 $(LIB)(dendro.o) \
 $(LIB)(denlst.o) \
 $(LIB)(denov.o)  \
 $(LIB)(deucl.o)  \
 $(LIB)(dhide.o)  \
 $(LIB)(diff1o.o) \
 $(LIB)(dilation.o) \
 $(LIB)(dirnm.o)   \
 $(LIB)(disc.o)    \
 $(LIB)(discrim.o) \
 $(LIB)(disp.o)    \
 $(LIB)(dist1.o)   \
 $(LIB)(dist2.o)   \
 $(LIB)(dist4.o)   \
 $(LIB)(dist_p.o)  \
 $(LIB)(distest.o) \
 $(LIB)(divcr3.o) \
 $(LIB)(dust.o) \
 $(LIB)(dvcr.o) \
 $(LIB)(docorr3_s.o) \
 $(LIB)(docors3_s.o) \
 $(LIB)(docseparate.o) \
 $(LIB)(docsub.o) \
 $(LIB)(docs1.o) \
 $(LIB)(domax3_s.o) \
 $(LIB)(domin3_s.o) \
 $(LIB)(dprofd.o) \
 $(LIB)(dprofl.o) \
 $(LIB)(driv1.o) \
 $(LIB)(driv2.o) \
 $(LIB)(driv3.o) \
 $(LIB)(edge.o) \
 $(LIB)(edm.o) \
 $(LIB)(eff.o) \
 $(LIB)(ehist.o) \
 $(LIB)(eigpercent.o) \
 $(LIB)(empsli.o) \
 $(LIB)(endit.o) \
 $(LIB)(ener.o) \
 $(LIB)(enfr_2.o) \
 $(LIB)(enhanc.o) \
 $(LIB)(entrop.o) \
 $(LIB)(envelope.o) \
 $(LIB)(epur4.o)   \
 $(LIB)(erc.o)     \
 $(LIB)(erosion.o) \
 $(LIB)(errf2.o) \
 $(LIB)(erri2.o) \
 $(LIB)(errt.o)  \
 $(LIB)(evalnq.o)\
 $(LIB)(evalsympar.o) \
 $(LIB)(express3q.o)  \
 $(LIB)(exprq.o)  \
 $(LIB)(falb.o)   \
 $(LIB)(fallg.o)  \
 $(LIB)(fbs3.o)   \
 $(LIB)(fbs2.o)   \
 $(LIB)(fcnq.o)   \
 $(LIB)(fermp.o)  \
 $(LIB)(fftc_d.o) \
 $(LIB)(fftc_q.o) \
 $(LIB)(fftmcf.o) \
 $(LIB)(fftr_d.o) \
 $(LIB)(fftr_q.o) \
 $(LIB)(fgr.o)    \
 $(LIB)(fhodt.o)  \
 $(LIB)(ffilts.o) \
 $(LIB)(fillbessi0.o) \
 $(LIB)(filcad.o) \
 $(LIB)(filcan.o) \
 $(LIB)(fildat.o) \
 $(LIB)(filelist.o) \
 $(LIB)(filerd.o) \
 $(LIB)(filgen.o) \
 $(LIB)(filget.o) \
 $(LIB)(filget_at.o) \
 $(LIB)(filnamandext.o) \
 $(LIB)(filnamsub.o) \
 $(LIB)(filseqp.o) \
 $(LIB)(filsli.o)  \
 $(LIB)(filtdoc.o) \
 $(LIB)(filter.o)  \
 $(LIB)(filter_har.o) \
 $(LIB)(filtplot.o) \
 $(LIB)(findfrq.o)  \
 $(LIB)(findcent.o) \
 $(LIB)(findfilcen.o) \
 $(LIB)(findlbq.o)  \
 $(LIB)(findminsormaxs.o) \
 $(LIB)(findmx.o)   \
 $(LIB)(findmx_f.o) \
 $(LIB)(findmx_q.o) \
 $(LIB)(findridges.o) \
 $(LIB)(fint.o)  \
 $(LIB)(fint3.o) \
 $(LIB)(fitt.o)  \
 $(LIB)(flushresults.o) \
 $(LIB)(flushfile.o) \
 $(LIB)(fmax_c.o) \
 $(LIB)(fmax_q.o) \
 $(LIB)(fmax_s.o) \
 $(LIB)(fmin_c.o) \
 $(LIB)(fmin_q.o) \
 $(LIB)(fmin_s.o) \
 $(LIB)(fmr_1.o)  \
 $(LIB)(fmr_2.o)  \
 $(LIB)(fmrs.o)   \
 $(LIB)(fmrs_1.o) \
 $(LIB)(fmrs_1d.o) \
 $(LIB)(fmrs_2.o)  \
 $(LIB)(fmrs_2dr.o) \
 $(LIB)(fmrs_2r.o) \
 $(LIB)(fmrs_3.o)  \
 $(LIB)(fmrs_3r.o) \
 $(LIB)(fmrs_3dr.o) \
 $(LIB)(fnrm.o) \
 $(LIB)(four1.o) \
 $(LIB)(four1a_fp.o) \
 $(LIB)(four1b.o) \
 $(LIB)(four1c.o) \
 $(LIB)(fouring.o) \
 $(LIB)(fouring_q.o) \
 $(LIB)(four_fq.o) \
 $(LIB)(frng.o) \
 $(LIB)(framealign_mlr.o) \
 $(LIB)(frngs.o) \
 $(LIB)(fromtoq.o) \
 $(LIB)(fsort.o) \
 $(LIB)(ftr.o) \
 $(LIB)(fv.o) \
 $(LIB)(gali.o) \
 $(LIB)(gam.o) \
 $(LIB)(gaussp.o) \
 $(LIB)(gee.o) \
 $(LIB)(getcoo.o) \
 $(LIB)(getcoot.o) \
 $(LIB)(getdefocus.o) \
 $(LIB)(getdocdat.o)  \
 $(LIB)(getdoclist.o) \
 $(LIB)(getfilenum.o) \
 $(LIB)(gethedccp4.o) \
 $(LIB)(getima_q.o) \
 $(LIB)(getlbno.o)  \
 $(LIB)(getlab.o)   \
 $(LIB)(get_tf_input.o) \
 $(LIB)(getnexttoken.o) \
 $(LIB)(getsym.o)   \
 $(LIB)(getthreads.o) \
 $(LIB)(gnc.o)   \
 $(LIB)(gnc2s.o) \
 $(LIB)(gotbranch.o) \
 $(LIB)(gp.o)   \
 $(LIB)(gprp.o) \
 $(LIB)(graphp.o) \
 $(LIB)(graphs.o) \
 $(LIB)(gsmod.o) \
 $(LIB)(gsmodl.o) \
 $(LIB)(hali.o) \
 $(LIB)(hclp.o) \
 $(LIB)(hcls.o) \
 $(LIB)(hdiag.o) \
 $(LIB)(hdls.o) \
 $(LIB)(hels.o) \
 $(LIB)(hgr.o) \
 $(LIB)(hiang.o) \
 $(LIB)(hismap.o) \
 $(LIB)(hismap4.o) \
 $(LIB)(hisd.o) \
 $(LIB)(hist.o) \
 $(LIB)(histe.o) \
 $(LIB)(histod.o) \
 $(LIB)(histog.o) \
 $(LIB)(histodc.o) \
 $(LIB)(hkmc.o) \
 $(LIB)(iftorpn.o) \
 $(LIB)(imseq4.o) \
 $(LIB)(imsq.o) \
 $(LIB)(imstat.o) \
 $(LIB)(imtovol.o) \
 $(LIB)(incor3.o) \
 $(LIB)(incort.o) \
 $(LIB)(inf3_p.o) \
 $(LIB)(inf_q.o) \
 $(LIB)(initunix.o) \
 $(LIB)(inlnbuf.o) \
 $(LIB)(inquirealloc.o) \
 $(LIB)(inquirearch.o) \
 $(LIB)(inquirecomp.o) \
 $(LIB)(inquirehead.o) \
 $(LIB)(inquireif.o) \
 $(LIB)(inquireif1.o) \
 $(LIB)(inquirereg.o) \
 $(LIB)(inquiresync.o) \
 $(LIB)(inttochar.o) \
 $(LIB)(interp.o) \
 $(LIB)(interps.o) \
 $(LIB)(interp3.o) \
 $(LIB)(interp_fbs3.o) \
 $(LIB)(interp_fbs.o) \
 $(LIB)(invbt.o) \
 $(LIB)(invgdist.o) \
 $(LIB)(ipalign64.o) \
 $(LIB)(irmrecord.o) \
 $(LIB)(irp3.o) \
 $(LIB)(ischar.o) \
 $(LIB)(isdigi.o) \
 $(LIB)(isbare.o) \
 $(LIB)(isdocincore.o) \
 $(LIB)(itpow.o) \
 $(LIB)(johist.o) \
 $(LIB)(johistf.o) \
 $(LIB)(jpmsk1.o) \
 $(LIB)(jpmsk2.o) \
 $(LIB)(jpmsk3.o) \
 $(LIB)(labelspi.o) \
 $(LIB)(lastoccur.o) \
 $(LIB)(latcen.o) \
 $(LIB)(latice.o) \
 $(LIB)(lenopenfile.o) \
 $(LIB)(listit.o) \
 $(LIB)(listitr.o) \
 $(LIB)(listhedmrc.o) \
 $(LIB)(lmdif.o) \
 $(LIB)(lmdif1.o) \
 $(LIB)(lnblnk.o) \
 $(LIB)(lnblnkn.o) \
 $(LIB)(local.o) \
 $(LIB)(lhist.o) \
 $(LIB)(log2.o) \
 $(LIB)(logifq.o) \
 $(LIB)(longer.o) \
 $(LIB)(macf.o) \
 $(LIB)(macf_p.o) \
 $(LIB)(mahal.o) \
 $(LIB)(make_close_list.o) \
 $(LIB)(maktab.o) \
 $(LIB)(malfi3.o) \
 $(LIB)(mapdist.o) \
 $(LIB)(mapdoc.o) \
 $(LIB)(mapfilt.o) \
 $(LIB)(mapim.o) \
 $(LIB)(mapstat.o) \
 $(LIB)(mask.o) \
 $(LIB)(maskmu.o) \
 $(LIB)(matinv.o) \
 $(LIB)(maxf.o) \
 $(LIB)(mccf.o) \
 $(LIB)(mccf_p.o) \
 $(LIB)(mcircl.o) \
 $(LIB)(md2.o) \
 $(LIB)(md3.o) \
 $(LIB)(meanshift.o) \
 $(LIB)(median.o) \
 $(LIB)(meed.o) \
 $(LIB)(mehist.o) \
 $(LIB)(mirror.o) \
 $(LIB)(mlc.o) \
 $(LIB)(mline.o) \
 $(LIB)(mline3.o) \
 $(LIB)(mod1.o) \
 $(LIB)(model.o) \
 $(LIB)(model3.o) \
 $(LIB)(mont.o) \
 $(LIB)(motif1.o) \
 $(LIB)(motif_init.o) \
 $(LIB)(motif_lsd.o) \
 $(LIB)(motif.o) \
 $(LIB)(mpi.o) \
 $(LIB)(mpoint.o) \
 $(LIB)(mr2to3d.o) \
 $(LIB)(mralign.o) \
 $(LIB)(mrand.o) \
 $(LIB)(mrang2.o) \
 $(LIB)(mrcalerr.o) \
 $(LIB)(mrcg3.o) \
 $(LIB)(mrcg3d.o) \
 $(LIB)(mrcp.o) \
 $(LIB)(mrdoc.o) \
 $(LIB)(mrdoccor.o) \
 $(LIB)(mread.o) \
 $(LIB)(mrectl.o) \
 $(LIB)(mreigen.o) \
 $(LIB)(mrerror.o) \
 $(LIB)(mrgetinfo.o) \
 $(LIB)(mrk.o) \
 $(LIB)(mrk3.o) \
 $(LIB)(mrkur3.o) \
 $(LIB)(mrncolor.o) \
 $(LIB)(mrnewangle.o) \
 $(LIB)(mrproj.o)   \
 $(LIB)(mrprrep.o)  \
 $(LIB)(mrputinfo.o)\
 $(LIB)(mrquater.o) \
 $(LIB)(mrrefl.o)   \
 $(LIB)(mrreplace.o)\
 $(LIB)(mrrotate.o) \
 $(LIB)(mrrsurf.o)  \
 $(LIB)(mrrt.o)     \
 $(LIB)(mrscale.o)  \
 $(LIB)(mrscale3.o)\
 $(LIB)(mrshift.o)\
 $(LIB)(mrsurf.o) \
 $(LIB)(mrtheta.o)\
 $(LIB)(msum.o)   \
 $(LIB)(mtnv.o)   \
 $(LIB)(mtpr.o)   \
 $(LIB)(mtrian.o) \
 $(LIB)(mult.o)   \
 $(LIB)(multishift.o) \
 $(LIB)(mx.o)       \
 $(LIB)(mygetenv.o) \
 $(LIB)(mymods.o) \
 $(LIB)(myread.o) \
 $(LIB)(mytime.o) \
 $(LIB)(mywr0.o)  \
 $(LIB)(mywr1.o)  \
 $(LIB)(mywr2.o)  \
 $(LIB)(narea.o)  \
 $(LIB)(negate.o) \
 $(LIB)(negati.o) \
 $(LIB)(newkmeans.o) \
 $(LIB)(nextversion.o) \
 $(LIB)(nn4.o) \
 $(LIB)(nn24.o) \
 $(LIB)(noise.o) \
 $(LIB)(norm3.o) \
 $(LIB)(normas.o) \
 $(LIB)(normasc.o) \
 $(LIB)(normm.o) \
 $(LIB)(norppf.o) \
 $(LIB)(noyau.o) \
 $(LIB)(nrml.o) \
 $(LIB)(numdig.o) \
 $(LIB)(od.o) \
 $(LIB)(oneline.o) \
 $(LIB)(opauxfile.o) \
 $(LIB)(openchk.o) \
 $(LIB)(opendoc.o) \
 $(LIB)(openfil.o) \
 $(LIB)(openfil_mrc.o) \
 $(LIB)(openfil_o_mrc.o) \
 $(LIB)(openfil_n_mrc.o) \
 $(LIB)(openinstk.o) \
 $(LIB)(openstk.o) \
 $(LIB)(opfilec.o) \
 $(LIB)(opfiles.o) \
 $(LIB)(opfiles_mrc.o) \
 $(LIB)(opstreamfile.o) \
 $(LIB)(opsf3_p.o) \
 $(LIB)(opsf_p.o) \
 $(LIB)(oracfmsk.o) \
 $(LIB)(oracfmskm.o) \
 $(LIB)(orcdf.o) \
 $(LIB)(orient.o) \
 $(LIB)(ormap.o) \
 $(LIB)(ormd.o) \
 $(LIB)(ormd_p.o) \
 $(LIB)(outim_q.o) \
 $(LIB)(outpr.o) \
 $(LIB)(outrng.o) \
 $(LIB)(parabl.o) \
 $(LIB)(parabld.o) \
 $(LIB)(parsesympar.o) \
 $(LIB)(parst.o) \
 $(LIB)(patch.o) \
 $(LIB)(pdates.o) \
 $(LIB)(pdb.o)     \
 $(LIB)(picksl.o)  \
 $(LIB)(picksv.o)  \
 $(LIB)(pj2d.o)    \
 $(LIB)(pj3g.o)    \
 $(LIB)(pj3_fbs.o) \
 $(LIB)(pj3_n.o)   \
 $(LIB)(pj3q_n.o)  \
 $(LIB)(pkd.o)     \
 $(LIB)(pksr3.o)   \
 $(LIB)(plnedg.o)  \
 $(LIB)(plot1.o)   \
 $(LIB)(plot2.o)   \
 $(LIB)(polar_cc.o) \
 $(LIB)(polish.o) \
 $(LIB)(poj.o)    \
 $(LIB)(polfit.o) \
 $(LIB)(pollsq.o) \
 $(LIB)(polqs.o) \
 $(LIB)(posaxis.o) \
 $(LIB)(postrt.o) \
 $(LIB)(proc_set.o) \
 $(LIB)(pr3d.o) \
 $(LIB)(pr3d_new.o) \
 $(LIB)(pr3db.o) \
 $(LIB)(prb1d.o) \
 $(LIB)(prb1.o) \
 $(LIB)(prepcub.o) \
 $(LIB)(prepcub_s.o) \
 $(LIB)(prepsl_2.o) \
 $(LIB)(prjc0.o) \
 $(LIB)(prjc2.o) \
 $(LIB)(prjcq.o) \
 $(LIB)(prjs2.o) \
 $(LIB)(prober.o) \
 $(LIB)(prp.o) \
 $(LIB)(prt.o) \
 $(LIB)(pttern.o) \
 $(LIB)(pts_on_sphere.o) \
 $(LIB)(putlin.o) \
 $(LIB)(putpt.o) \
 $(LIB)(putpt1.o) \
 $(LIB)(putpt2.o) \
 $(LIB)(putyslice.o) \
 $(LIB)(putzslice.o) \
 $(LIB)(pw2sdr.o) \
 $(LIB)(pw2sr.o) \
 $(LIB)(pw3sdr.o) \
 $(LIB)(pw3sr.o) \
 $(LIB)(qali.o) \
 $(LIB)(qfact.o) \
 $(LIB)(qfact_p.o) \
 $(LIB)(qnrf.o) \
 $(LIB)(qphi.o) \
 $(LIB)(qstat.o) \
 $(LIB)(qstatloc.o) \
 $(LIB)(qtheta.o) \
 $(LIB)(quadri.o) \
 $(LIB)(quadri_q.o) \
 $(LIB)(radav.o) \
 $(LIB)(radavi.o) \
 $(LIB)(ramp_p.o) \
 $(LIB)(rand_p.o) \
 $(LIB)(rann.o) \
 $(LIB)(rathin.o) \
 $(LIB)(rawtospider.o) \
 $(LIB)(rclu.o) \
 $(LIB)(rcnv2_p.o) \
 $(LIB)(rcnv3_p.o) \
 $(LIB)(rconv.o) \
 $(LIB)(rctfss.o) \
 $(LIB)(rctfone.o) \
 $(LIB)(rdclus.o) \
 $(LIB)(rddocn.o) \
 $(LIB)(udoc.o) \
 $(LIB)(rdpa.o) \
 $(LIB)(rdpr.o) \
 $(LIB)(rdpra.o) \
 $(LIB)(rdprai.o) \
 $(LIB)(rdpranc.o) \
 $(LIB)(rdprinc.o) \
 $(LIB)(rdpris.o) \
 $(LIB)(rdpri1s.o) \
 $(LIB)(rdpri2s.o) \
 $(LIB)(rdpri3s.o) \
 $(LIB)(rdprm.o)   \
 $(LIB)(rdprm2.o)  \
 $(LIB)(rdprm1s.o) \
 $(LIB)(rdprm2s.o) \
 $(LIB)(rdprm3s.o) \
 $(LIB)(rdprmc.o) \
 $(LIB)(rdprmi.o) \
 $(LIB)(rdprop.o) \
 $(LIB)(rdsl_p.o) \
 $(LIB)(readpdb.o) \
 $(LIB)(readrq.o) \
 $(LIB)(reads.o) \
 $(LIB)(readv.o) \
 $(LIB)(regparse.o) \
 $(LIB)(reg_set.o) \
 $(LIB)(redlin.o) \
 $(LIB)(redlin8.o) \
 $(LIB)(rednpadvol.o) \
 $(LIB)(redpro2.o) \
 $(LIB)(redrec8.o) \
 $(LIB)(redvol.o) \
 $(LIB)(reform.o) \
 $(LIB)(reform0.o) \
 $(LIB)(remez.o) \
 $(LIB)(remezp.o) \
 $(LIB)(bpcg.o) \
 $(LIB)(repr2_s.o) \
 $(LIB)(retir.o) \
 $(LIB)(reversebytes.o) \
 $(LIB)(rew.o) \
 $(LIB)(rewf.o) \
 $(LIB)(rfactsd2.o) \
 $(LIB)(rfactsdo.o) \
 $(LIB)(rgri.o) \
 $(LIB)(ringwe.o) \
 $(LIB)(rmfilterz.o) \
 $(LIB)(rmfilty.o) \
 $(LIB)(rmfiltz.o) \
 $(LIB)(rmiwf.o) \
 $(LIB)(rmmakemask.o) \
 $(LIB)(rmrepl.o)     \
 $(LIB)(rmz_p.o)      \
 $(LIB)(rotates.o)    \
 $(LIB)(rotates3l.o)  \
 $(LIB)(rotates3q.o)  \
 $(LIB)(rotatesl3q.o) \
 $(LIB)(rotates3fbs.o) \
 $(LIB)(rot32.o) \
 $(LIB)(rot3_p.o) \
 $(LIB)(rot_p.o) \
 $(LIB)(rotal3.o) \
 $(LIB)(rotas3.o) \
 $(LIB)(rotl3.o) \
 $(LIB)(rotl3q.o) \
 $(LIB)(rotqss.o) \
 $(LIB)(rots3.o) \
 $(LIB)(rots3q.o) \
 $(LIB)(rprq.o) \
 $(LIB)(rtsf.o) \
 $(LIB)(rtsq.o) \
 $(LIB)(rtq.o) \
 $(LIB)(rtq_q.o) \
 $(LIB)(rtqs_q.o) \
 $(LIB)(rys1.o) \
 $(LIB)(rys3.o) \
 $(LIB)(saqb.o) \
 $(LIB)(savd.o) \
 $(LIB)(savdn1.o) \
 $(LIB)(savdocq.o) \
 $(LIB)(sclassi.o) \
 $(LIB)(scoran3.o) \
 $(LIB)(searchq.o) \
 $(LIB)(secnds.o)  \
 $(LIB)(seedfill.o) \
 $(LIB)(seeds.o) \
 $(LIB)(sen3a.o) \
 $(LIB)(setenv.o) \
 $(LIB)(sethedccp4.o) \
 $(LIB)(sethead.o) \
 $(LIB)(setlab.o)  \
 $(LIB)(setmode.o) \
 $(LIB)(setmode_mrc.o) \
 $(LIB)(setnumrings.o) \
 $(LIB)(setpart.o) \
 $(LIB)(setprm.o)  \
 $(LIB)(setprmb.o) \
 $(LIB)(setsympar.o) \
 $(LIB)(setthreads.o) \
 $(LIB)(setval.o) \
 $(LIB)(sdocic.o) \
 $(LIB)(sgraf.o)  \
 $(LIB)(sh180_2.o) \
 $(LIB)(shelk.o)  \
 $(LIB)(shfc_2.o) \
 $(LIB)(shfi_2.o) \
 $(LIB)(shfi_q.o) \
 $(LIB)(shfm_2.o) \
 $(LIB)(shift1.o) \
 $(LIB)(shift2.o) \
 $(LIB)(shift3.o) \
 $(LIB)(shift_3d.o) \
 $(LIB)(shift_pf.o) \
 $(LIB)(shiftr.o) \
 $(LIB)(shifts.o) \
 $(LIB)(shiftt.o) \
 $(LIB)(shosli.o) \
 $(LIB)(shosli2.o) \
 $(LIB)(shrinkq.o) \
 $(LIB)(shuffledoc.o) \
 $(LIB)(single.o) \
 $(LIB)(sizchk.o) \
 $(LIB)(skeleton.o) \
 $(LIB)(slice.o)  \
 $(LIB)(smt3.o)   \
 $(LIB)(smt3_q.o) \
 $(LIB)(snrb.o)   \
 $(LIB)(snrf.o)   \
 $(LIB)(softmask.o) \
 $(LIB)(solv2d.o) \
 $(LIB)(solve.o)  \
 $(LIB)(sort.o)   \
 $(LIB)(sorti.o)  \
 $(LIB)(sortint.o)\
 $(LIB)(sortz.o)  \
 $(LIB)(spdtotiff.o) \
 $(LIB)(sphdecon.o) \
 $(LIB)(spca3.o)  \
 $(LIB)(speak.o)  \
 $(LIB)(speak3.o) \
 $(LIB)(speakc.o) \
 $(LIB)(speakm.o) \
 $(LIB)(speakr.o) \
 $(LIB)(spireout.o) \
 $(LIB)(spotwt.o) \
 $(LIB)(squeez.o) \
 $(LIB)(ssnrb.o)  \
 $(LIB)(ssnr3.o)  \
 $(LIB)(ssnr3dnn.o) \
 $(LIB)(sspush.o) \
 $(LIB)(sslowcas.o) \
 $(LIB)(ssupcas.o) \
 $(LIB)(stabk.o) \
 $(LIB)(stack.o) \
 $(LIB)(star2doc.o) \
 $(LIB)(step.o) \
 $(LIB)(stepi.o) \
 $(LIB)(subaf.o) \
 $(LIB)(subchar.o) \
 $(LIB)(subkmeans.o) \
 $(LIB)(subkmns.o) \
 $(LIB)(subsympar.o) \
 $(LIB)(sumali.o) \
 $(LIB)(sumali3.o) \
 $(LIB)(sumap.o)     \
 $(LIB)(sumeuler.o)  \
 $(LIB)(sunonly.o)   \
 $(LIB)(surfarea.o)  \
 $(LIB)(surffit.o)   \
 $(LIB)(surfcomp.o)  \
 $(LIB)(surftovol.o) \
 $(LIB)(symang.o)  \
 $(LIB)(sympar.o)  \
 $(LIB)(symvol.o)  \
 $(LIB)(tcnp.o)    \
 $(LIB)(tf.o)      \
 $(LIB)(tfcrf.o)   \
 $(LIB)(tfd.o)     \
 $(LIB)(tfed.o)    \
 $(LIB)(tflm4.o)   \
 $(LIB)(tfsnr.o)   \
 $(LIB)(thresh.o)  \
 $(LIB)(tilt.o)    \
 $(LIB)(tiltfd.o)  \
 $(LIB)(tima.o)    \
 $(LIB)(title.o)   \
 $(LIB)(toma.o)    \
 $(LIB)(to_peaks.o) \
 $(LIB)(to_polar.o) \
 $(LIB)(to_rays.o) \
 $(LIB)(traf.o)    \
 $(LIB)(trafc.o)   \
 $(LIB)(trafc3.o)  \
 $(LIB)(trafd.o)   \
 $(LIB)(trafl.o)   \
 $(LIB)(trafsim.o) \
 $(LIB)(transpose.o) \
 $(LIB)(tridi.o)   \
 $(LIB)(trinter.o) \
 $(LIB)(triquad.o) \
 $(LIB)(tswitch.o) \
 $(LIB)(ttest.o)   \
 $(LIB)(ttest1.o)  \
 $(LIB)(tvwn3.o)   \
 $(LIB)(udocic.o)  \
 $(LIB)(unixtounix8.o) \
 $(LIB)(unsav.o)    \
 $(LIB)(unsdal.o)   \
 $(LIB)(updf.o)     \
 $(LIB)(updtc.o)    \
 $(LIB)(updtf.o)    \
 $(LIB)(updtf_r.o)  \
 $(LIB)(uqu.o)      \
 $(LIB)(util1.o)    \
 $(LIB)(util2.o)    \
 $(LIB)(util2sup.o) \
 $(LIB)(util3.o)    \
 $(LIB)(util4.o)    \
 $(LIB)(util5.o)    \
 $(LIB)(util6.o)    \
 $(LIB)(util7.o)    \
 $(LIB)(util_11.o)   \
 $(LIB)(util_1010.o) \
 $(LIB)(util_1011.o) \
 $(LIB)(util_1110.o) \
 $(LIB)(varf.o)      \
 $(LIB)(varf1.o)     \
 $(LIB)(var3r.o)     \
 $(LIB)(var3d.o)     \
 $(LIB)(vax32u.o)    \
 $(LIB)(vaxtounix.o) \
 $(LIB)(vismap.o) \
 $(LIB)(vms.o)    \
 $(LIB)(vms_cd.o) \
 $(LIB)(voda.o)   \
 $(LIB)(voea.o)   \
 $(LIB)(vonea.o)  \
 $(LIB)(voept.o)  \
 $(LIB)(voia.o)   \
 $(LIB)(vomd.o)   \
 $(LIB)(vomq.o)   \
 $(LIB)(vora.o)   \
 $(LIB)(vprop.o)  \
 $(LIB)(vtil2.o)  \
 $(LIB)(vtil3.o)  \
 $(LIB)(warp3.o)  \
 $(LIB)(wate.o)   \
 $(LIB)(watershed.o) \
 $(LIB)(watershd.o)  \
 $(LIB)(wftcirc.o) \
 $(LIB)(wgbp2.o)   \
 $(LIB)(winave.o)  \
 $(LIB)(winave2.o) \
 $(LIB)(window.o)  \
 $(LIB)(bp3f.o)    \
 $(LIB)(wiw3g.o)   \
 $(LIB)(bp32f.o)   \
 $(LIB)(wpdp.o)      \
 $(LIB)(wpro_n.o)    \
 $(LIB)(wpro_fbs.o)  \
 $(LIB)(writepict.o) \
 $(LIB)(writepictn.o)\
 $(LIB)(writev.o)    \
 $(LIB)(writpro_n.o) \
 $(LIB)(wrsl_p.o)    \
 $(LIB)(wrtlin.o)    \
 $(LIB)(wrtlin8.o)   \
 $(LIB)(wrtseq8.o)   \
 $(LIB)(wrtvol.o)    \
 $(LIB)(wtf.o)       \
 $(LIB)(wtm.o)       \
 $(LIB)(xmipp.o)     \
 $(LIB)(xmselfile.o) 


